// -*-c++-*-
#ifndef STORMM_LAYERED_POTENTIAL_H
#define STORMM_LAYERED_POTENTIAL_H

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Math/formulas.h"
#include "Math/log_scale_spline.h"
#include "Math/matrix_ops.h"
#include "Potential/pme_util.h"
#include "Structure/structure_enumerators.h"
#include "energy_enumerators.h"
#include "layered_potential_metrics.h"

namespace stormm {
namespace energy {

using card::Hybrid;
using card::HybridTargetLevel;
using energy::elecPMEDirectSpace;
using stmath::sigmoid;
using stmath::sigmoidf;
using stmath::LogScaleSpline;
using stmath::LogSplineTable;
using stmath::qrSolver;

/// \brief The abstract of the LayeredPotential class provides only one of the layers.  Kernels
///        can be called with one or more such abstracts, but due to the logarithmic nature of
///        the volume of particles applying each potential a "one size fits all" implementation
///        is not practical.  Furthermore, due to the way in which successive layers of the
///        potential soften, multiple time step approaches which only invoke one to two layers of
///        the potentials at a time are feasible.
template <typename T4>
class NrgLayerKit {

  /// The constructor takes inputs for all parameters, the bulk of which is the LogScaleSpline
  /// abstract.  
  NrgLayerKit(const LogSplineTable<T4> &u_in, const LogSplineTable<T4> &du_in,
              const LogSplineTable<T4> &d2u_in, const LogSplineTable<T4> &d3u_in,
              double cutoff_in);

  /// \brief The default copy and move constructors are valid, but the presence of const members
  ///        (which even contain other const members, in the case of the LogSplineTable abstracts)
  ///        prohibits copy and move assignment operators.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object placed on the right hand side of the assignment statement
  /// \{
  NrgLayerKit(const NrgLayerKit &original) = default;
  NrgLayerKit(NrgLayerKit &&original) = default;
  NrgLayerKit& operator=(const NrgLayerKit &original) = default;
  NrgLayerKit& operator=(NrgLayerKit &&original) = default;
  /// \}

  // A total of four logarithmic spline tables are contained within the abstract.
  const LogSplineTable<T4> u;    ///< The potential
  const LogSplineTable<T4> du;   ///< The first derivative
  const LogSplineTable<T4> d2u;  ///< The second derivative
  const LogSplineTable<T4> d3u;  ///< The third derivative
  const double cutoff;           ///< The cutoff at which the layer stops completely.  It may be
                                 ///<   advantageous to apply all interactions regardless of the
                                 ///<   cutoff, as the spline tables will register zero for all
                                 ///<   segment beyond the cutoff.  However, spline segments that
                                 ///<   include the cutoff will only go to approximately zero at
                                 ///<   the cutoff, a fact which may be particularly relevant for
                                 ///<   the third derivative which can approach zero with a
                                 ///<   considerable slope at the cutoff.
};
  
/// \brief A long-ranged potential can be broken down into a series of successively smoother
///        potentials applicable over longer and longer ranges.  This object will manage a process
///        whereby each potential and its first three derivatives decay exactly to zero at the
///        specified boundary.  The process generates some very complicated forms, particularly in
///        the shortest-range potential where the effects assigned to all longer-ranged potentials
///        must be subtracted off, but use of logarithmic adaptive splines can absorb the
///        complexity.
template <typename T, typename T4> class LayeredPotential {
public:

  /// \brief The constructor requires an enumeration to specify a particular potential form as
  ///        well as a rough understanding of the box size.  The CLASH form corresponds to an
  ///        empty object.  For van-der Waals potentials, the type of mixing rule can also be
  ///        significant: the shortest-ranged potentials will turn on the 1/r^6 potential with a
  ///        sigmoidal switching function to allow NBFIX and Lorentz-Berthelot style mixing rules
  ///        to apply their unique interactions at close range, pushing back the transition to a
  ///        geometric combining rule to a longer range where the difference between mixing rules
  ///        is trivial.
  /// \{
  LayeredPotential(const LayeredPotentialMetrics &parameters_in = LayeredPotentialMetrics());
  /// \}

  /// \brief So long as the nested parameters manager and LogScaleSpline objects copy and move
  ///        correctly, the copy and move constructors, as well as the copy and move assignment
  ///        operators will be valid.
  ///
  /// \param original  The existing object to copy or move
  /// \param other     Another object to place on the right hand side of the assignment statement
  /// \{
  LayeredPotential(const LayeredPotential &original) = default;
  LayeredPotential(LayeredPotential &&original) = default;
  LayeredPotential& operator=(const LayeredPotential &other) = default;
  LayeredPotential& operator=(LayeredPotential &&other) = default;
  /// \}

  /// \brief Get a const reference to the underlying metrics, to view specifications such as the
  ///        cutoff at a particular level, one of the exponential scaling factors, etc.  All
  ///        accessors in the nested class can then be utilized without copying the code.
  const LayeredPotentialMetrics& getParameters() const;

  /// \brief Obtain the quartet of coefficients for one of the layers.
  ///
  /// \param layer  Index of the layer of interest
  double4 getSmoothingCoefficients(int layer) const;
  
  /// \brief Get the value of one of the potential functions at a specified layer and distance,
  ///        using the appropriate internal spline table.
  ///
  /// \param layer  Index of the layer of interest
  /// \param r      Distance at which to evaluate the interaction (specifying a value beyond the
  ///               particular layer's cutoff will return zero)
  T getSplinedValue(int layer, T r, T r2 = -1.0) const;

  /// \brief Compute the value of one of the potential functions, using analytic expressions, at a
  ///        specified layer and distance.  Descriptions of input parameters follow from
  ///        getSplinedValue() above.
  T getAnalyticValue(int layer, T r, T r2 = -1.0) const;

  /// \brief Get the value of the potential derivative at a specified layer and distance, making
  ///        use of an internal table for the derivative.
  /// 
  /// \param layer  Index of the layer of interest
  /// \param r      Distance at which to evaluate the interaction (specifying a value beyond the
  ///               particular layer's cutoff will return zero)
  /// \param order  The order of the derivative (default 1, the first derivative--specifying zero
  ///               will return the function value)
  T getSplinedDerivative(int layer, T r, T r2 = -1.0, int order = 1) const;

  /// \brief Get the value of the potential derivative at a specified layer and distance, using
  ///        analytic calculations.  Descriptions of input parameters follow from
  ///        getSplinedDerivative() above.
  T getAnalyticDerivative(int layer, T r, T r2 = -1.0, int order = 1) const;

  /// \brief Get a read-only abstract for the object, gear towards a specific layer of the
  ///        potential.  The abstract for the appropriate spline table, plus limits of the layer's
  ///        range, will be included.  Because the final spline segment might not terminate at
  ///        precisely the layer's cutoff, it is essential to check whether any particular
  ///        interaction is in the applicable range, rather than relying on spline table entries
  ///        to be zero past the applicable range.
  ///
  /// \param layer_index  Index of the layer for which to prepare the abstract
  /// \param tier         Assign pointers to target memory on the CPU host or GPU device
  const NrgLayerKit<T4> data(int layer_index,
                             HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  
private:
  LayeredPotentialMetrics parameters;
  std::vector<LogScaleSpline<T4>> spline_tables;
  std::vector<double4> smoothing_coefficients;

  /// \brief Check that a requested layer index is within the number of expected layers.
  ///
  /// \param layer_index_in  Index of the layer of interest
  /// \param caller          Name of the calling function, for error tracing
  void validateLayerIndex(int layer_index, const char* caller = nullptr) const;
  
  /// \brief Solve the exponential series coefficients for the decomposition of long-ranged
  ///        interactions.  See documentation in the function for more details of each level's
  ///        coefficients.
  void solveLayerCoefficients();
};

} // namespace energy
} // namespace stormm

#include "layered_potential.tpp"

#endif
