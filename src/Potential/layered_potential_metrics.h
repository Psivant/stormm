// -*-c++-*-
#ifndef STORMM_LAYERED_POTENTIAL_METRICS_H
#define STORMM_LAYERED_POTENTIAL_METRICS_H

#include <vector>
#include "copyright.h"
#include "Constants/scaling.h"
#include "Constants/symbol_values.h"
#include "Math/log_scale_spline.h"
#include "Structure/structure_enumerators.h"
#include "energy_enumerators.h"

namespace stormm {
namespace energy {

using structure::BoundaryCondition;
using symbols::amber_ancient_bioq;
  
/// \brief Constants to put bounds on the LayeredPotential
/// \{
constexpr double hail_minimum_range_limit = 4.0;
constexpr double hail_minimum_range_multiplier = 1.5;
constexpr double hail_maximum_range_multiplier = 3.0;
constexpr double hail_default_vdw_xtn_midpoint = 7.0;
constexpr double hail_default_vdw_xtn_intensity = 3.0;
constexpr double hail_default_ewald_coefficient = 0.3;
/// \}
  
/// \brief Collect the details needed to formulate a layered potential, including the length of
///        each layer's cutoff (for direct particle-particle interactions), expected size of the
///        unit cell (if applicable), mixing rule or Coulomb's constant, and the primary
///        exponential coefficients for fitting successive smoothing functions.
class LayeredPotentialMetrics {
public:

  /// \brief The constructor takes input for all member variables, but various modifier functions
  ///        can set them one at a time for clarity.
  ///
  /// \param particle_particle_cutoff  The cutoff that will be applied to interactions at the
  ///                                  lowest level of the potential, where there is no mesh and
  ///                                  all calculations are done with pairwise calculations between
  ///                                  particles
  LayeredPotentialMetrics(DecomposablePotential form_in = DecomposablePotential::ELECTROSTATIC,
                          BoundaryCondition boundaries_in = BoundaryCondition::ISOLATED,
                          double coulomb_in = amber_ancient_bioq,
                          VdwCombiningRule mixing_rule_in = VdwCombiningRule::GEOMETRIC,
                          const std::vector<double> &exponent_c_in = { 1.0, 1.5, 2.0 },
                          double particle_particle_cutoff = 5.0,
                          double range_multiplier_in = 2.0,
                          double range_limit_in = 1000.0,
                          double vdw_transition_midpoint_in = hail_default_vdw_xtn_midpoint,
                          double vdw_transition_intensity_in = hail_default_vdw_xtn_intensity,
                          double ewald_coefficient_in = hail_default_ewald_coefficient,
                          const std::vector<double> &box_vectors_in = { 50.0,  0.0,  0.0,
                                                                         0.0, 50.0,  0.0,
                                                                         0.0,  0.0, 50.0 });

  /// \brief Using Standard Template Library objects for anything that is not a scalar value, and
  ///        with no const members, the default copy and move constructors as well as copy and
  ///        move assignment operators apply.
  /// \{
  LayeredPotentialMetrics(const LayeredPotentialMetrics &original) = default;
  LayeredPotentialMetrics(LayeredPotentialMetrics &&original) = default;
  LayeredPotentialMetrics& operator=(const LayeredPotentialMetrics &other) = default;
  LayeredPotentialMetrics& operator=(LayeredPotentialMetrics &&other) = default;
  /// \}
  
  /// \brief Get the requested potential form.
  DecomposablePotential getForm() const;

  /// \brief Get the expected boundary conditions.
  BoundaryCondition getBoundaryCondition() const;
  
  /// \brief Get the number of layers by which the potential is split.
  int getLayerCount() const;
  
  /// \brief Get the cutoff value for one of the layers.  The cutoff for one layer is the handoff
  ///        point for the next layer up.
  ///
  /// \param layer_index  Index of the layer of interest
  double getCutoff(int layer_index) const;

  /// \brief Get the value of Coulomb's constant to use in preparing potentials.  This will raise
  ///        an exception if the form of the potential is not electrostatic.
  double getCoulombConstant() const;

  /// \brief Get the mixing rule applicable to dispersion interactions.  Specific Lennard-Jones
  ///        parameters are not needed at this stage--they will play a role when constructing the
  ///        LayeredPotential object based on the parameters in this preparatory class.
  VdwCombiningRule getMixingRule() const;

  /// \brief Get one of the three baseline exponential coefficients.  The strategy with Layered
  ///        Interpolation in Real Space is to split the potential at a given range Rh using a
  ///        softening function of the form:
  ///
  ///        U(r) = A exp(c1 * (r - Rh)) + B exp(c2 * (r - Rh)) + C exp(c3 * (r - Rh)) + D
  ///
  ///        Here, Rh is the "handoff" distance and is equal to the cutoff of the next lower layer
  ///        of the interaction, the lowest layer having the shortest cutoff and all interactions
  ///        being direct particle-particle pairs.  Potentials at all levels go to zero exactly,
  ///        with zero derivative, second derivative, and third derivative, at their own cutoff
  ///        range which coincides with handoff point Rh for the next layer.  The exponential
  ///        coefficients c1, c2, and c3 are parameters of the model which tune the rate at which
  ///        the smoothing potential moderates the original form.  {1.0, 1.5, 2.0} or
  ///        {1.25, 1.75, 2.25} are good choices for the electrostatic function and a handoff point
  ///        of about 5.0 Angstroms.
  ///
  ///        At successive layers of the potential, the exponential scaling factors are scaled down
  ///        as a proportion of the handoff point.  This does more than merely allow a single set
  ///        of coefficients to control the entire series.  It also helps maintain a consistent
  ///        stencil for mapping the finite potentials to meshes at a given level of accuracy.
  ///
  /// \param factor_index  Specify the index for c1, c2, or c3.  Values outside [0, 2] raise
  ///                      exceptions.
  /// \param layer_index   The layer of interest.  The lowet layer (index 0) is for direct
  ///                      particle-particle interactions and has no exponential series.
  ///                      Requesting such parameters will return 0.0.
  double getExponentFactor(int factor_index, int layer_index = 1) const;

  /// \brief Get the general scaling factor for the handoff points between successive layers.
  ///        Two is a natural choice, doubling the effective range of the potential in each layer,
  ///        but intermediate values such as 1.5 are also possible.  This, the particle-particle
  ///        interaction cutoff, and the stencil size (which is set for the simulation, outside of
  ///        the LayeredPotential) are the primary tunable parameters determining the accuracy of
  ///        the approximation.  The exponential factors are more minor parameters, in the sense
  ///        that they can make the approximation worse but there are generally "best" settings
  ///        given the choice of the particle-particle cutoff.
  double getRangeCompounding() const;

  /// \brief Get the maximum range for calculated interactions, applicable to systems with isolated
  ///        boundary conditions.
  double getMaximumRange() const;

  /// \brief Get the dispersion interaction switching midpoint, at which point whatever applicable
  ///        interaction is mixed 50:50 with the interaction expected for a geometric combining
  ///        rule.
  double getVdwTransitionMidpoint() const;

  /// \brief Get the intensity with which the dispersion interaction switches between its native
  ///        form and the interaction expected for a geometric combining rule.  The switching
  ///        function turns off the native potential by S(r) = 1 / (exp(p (r - r0)) + 1), where p
  ///        is the intensity returned by this function and r0 is the midpoint returned by
  ///        getVdwTransitionMidpoint() above.
  double getVdwTransitionIntensity() const;

  /// \brief Get the Ewald coefficient for decomposition of a PME direct space interaction.
  double getEwaldCoefficient() const;
  
  /// \brief Get the rough unit cell dimensions, expressed as a 3 x 3 matrix whose columns are the
  ///        Cartesian X, Y, and Z coordinates of each unit cell vector (it is identical to the
  ///        inverse transformation matrix for taking fractional coordinates into real space).
  const std::vector<double>& getUnitCellVectors() const;

  /// \brief Set the form of the potential that will be expressed.
  ///
  /// \param form_in  The potential to subdivide into layers
  void setForm(DecomposablePotential form_in);

  /// \brief Set the boundary conditions to use.
  ///
  /// \param boundary_in
  void setBoundaryCondition(BoundaryCondition boundary_in);

  /// \brief Set the cutoff for particle-particle interactions.
  ///
  /// \param cutoff_in  The chosen cutoff
  void setCutoff(double cutoff_in);

  /// \brief Set the value of Coulomb's constant.
  ///
  /// \param coulomb_in  The chosen value
  void setCoulombConstant(double coulomb_in);

  /// \brief Set the mixing rule to apply to Lennard-Jones parameter pairs.  This will determine
  ///        the forms of particle-particle interactions and the first grid layer potential.
  ///
  /// \param mixing_rule_in  The Lennard-Jones combination rule to apply
  void setMixingRule(VdwCombiningRule mixing_rule_in);

  /// \brief Set one of the baseline exponential factors.  These apply to the first smoothed
  ///        potential, the first layer involving a mesh after the particle-particle interactions.
  ///        Subsequent layers' exponential factors are scaled with each cutoff.
  ///
  /// Overloaded:
  ///   - Provide exponent factors one at a time, with term indices in each case
  ///   - Provide three exponent factors as a standard Template Library vector
  ///   - Provide three exponent factors as a double-precision tuple
  ///
  /// \param factor_index  Specify the index for c1, c2, or c3.  Values outside [0, 2] raise
  ///                      exceptions.
  /// \param exp_c_in      The exponential multiplier
  /// \{
  void setExponentFactor(int factor_index, double exp_c_in);
  void setExponentFactor(const std::vector<double> &exp_c_in);
  void setExponentFactor(const double3 exp_c_in);
  /// \}

  /// \brief Set the range compounding factor.
  ///
  /// \param range_multiplier_in  The compounding factor determining how much the cutoff increases
  ///                             from one layer of the approximation to the next
  void setRangeCompounding(double range_multiplier_in);

  /// \brief Set the maximum range for the layered potential cutoffs.  This will provide an upper
  ///        bound on the interactions in the case of isolated boundary conditions.
  ///
  /// \param range_limit_in  The maximum range for interactions by any of the smoothed potentials
  void setMaximumRange(double range_limit_in);
  
  /// \brief Set the dispersion interaction switching midpoint, at which point whatever applicable
  ///        interaction is mixed 50:50 with the interaction expected for a geometric combining
  ///        rule.
  ///
  /// \param midpoint_in  The transition midpoint to set
  void setVdwTransitionMidpoint(double midpoint_in);

  /// \brief Set the intensity with which the dispersion interaction switches between its native
  ///        form and the interaction expected for a geometric combining rule.  See the associated
  ///        accessor for details of the transition functional form..
  ///
  /// \param intensity_in  The transition intensity to set
  void setVdwTransitionIntensity(double intensity_in);

  /// \brief Set the Ewald coefficient for decomposition of a PME direct space interaction.
  ///
  /// \param ewald_coefficient_in  The chose Ewald coefficient, equal to 1/(2g) where g is the
  ///                              width of the Gaussian density smoothing form in the Ewald
  ///                              approximation
  void setEwaldCoefficient(double ewald_coefficient_in);
  
  /// \brief Set the unit cell vectors.  This will provide an upper bound on the interactions in
  ///        the case of periodic boundary conditions.
  ///
  /// \param box_vectors_in  The inverse transformation matrix taking particles from unit cell
  ///                        fractional space back to real space
  void setUnitCellVectors(const std::vector<double> &box_vectors_in);

private:

  DecomposablePotential form;       ///< Form of the potential being subdivided into layers
  BoundaryCondition boundaries;     ///< The type of boundary conditions that the potential should
                                    ///<   be prepared for.  A number of layers will be prepared
                                    ///<   that is sufficient to span the simulation box and, if
                                    ///<   necessary, include the effects of infinite lattice
                                    ///<   images.
  int layer_count;                  ///< The number of layers into which the potential subdivides
  std::vector<double> cutoffs;      ///< Cutoff points for each layer of the potential.  The
                                    ///<   cutoff of layer i is the handoff point at which layer
                                    ///<   i + 1 transitions to a softcore potential based on a
                                    ///<   three-term exponential series.
  double coulomb_constant;          ///< The value of Coulomb's constant defining an electrostatic
                                    ///<   potential (only valid for electrostatics)
  VdwCombiningRule mixing_rule;     ///< The manner in which different Lennard-Jones parameters
                                    ///<   for an atom type pair combine to define the pair
                                    ///<   interaction.  This affects how the lowest level of the
                                    ///<   approximation handles particle-particle interactions.
  std::vector<double3> exponent_c;  ///< Baseline exponential factors for three terms in the
                                    ///<   exponential series constituting the smoothing function
                                    ///<   for the first mesh-based potential (layer index 1).
                                    ///<   Coefficients for exponential terms in higher layers k
                                    ///<   are obtained by multiplying the values of exponent_c by
                                    ///<   cutoffs(1)/cutoffs(k).
  double range_multiplier;          ///< Range multipler for the cutoffs of successive levels.
                                    ///<   the final cutoff may break this series so as to become
                                    ///<   half the approximate width of the box along its shortest
                                    ///<   exit.
  double range_limit;               ///< The maximum range for any of the smoothed potentials.
                                    ///<   Interactions occurring at ranges longer than this could
                                    ///<   be missed if applying the layered potential in isolated
                                    ///<   boundary conditions.
  double vdw_transition_midpoint;   ///< The midpoint at which dispersion interactions will be
                                    ///<   mixed 50:50 with the interaction that would result from
                                    ///<   a geometric combination of the underlying parameters
  double vdw_transition_intensity;  ///< Controls the rate at which the long-ranged van-der Waals
                                    ///<   dispersion interaction transitions into a model based
                                    ///<   on the geometric parameter combining rule
  double ewald_coefficient;         ///< The Ewald coefficient to apply, in case the potential is
                                    ///<   decomposing some PME direct space interaction
  std::vector<double> box_vectors;  ///< Matrix of three column vectors defining the approximate
                                    ///<   size of the unit cell (for periodic as well as
                                    ///<   non-periodic boundary conditions--

  /// \brief Check that a requested layer index is within the number of expected layers.
  ///
  /// \param layer_index_in  Index of the layer of interest
  /// \param caller          Name of the calling function, for error tracing
  void validateLayerIndex(int layer_index, const char* caller = nullptr) const;

  /// \brief Compute the layer boundaries given the length of the cutoff on particle-particle
  ///        interactions, some notion of the range to which pair interactions will be necessary,
  ///        i.e. known boundary conditions plus the corresponding simulation size, and the
  ///        multiplicative factor between the cutoffs at successive levels.
  void computeLayers();
};

} // namespace energy
} // namespace stormm

#endif
