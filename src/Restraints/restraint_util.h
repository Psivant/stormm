// -*-c++-*-
#ifndef STORMM_RESTRAINT_UTIL_H
#define STORMM_RESTRAINT_UTIL_H

#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"

namespace stormm {
namespace restraints {

/// \brief Compute critical elements of the restraining potential: its difference from the target
///        value that determines some harmonic stiffness penalty, the harmonic penalty stiffness,
///        and the energy contribution.
///
/// \param init_k   Initial stiffness parameters
/// \param final_k  Final stiffness parameters
/// \param init_r   Initial displacement parameters
/// \param final_r  Final displacement parameters
/// \param mixwt    Pre-calculated mixing factor for combining initial and final parameters
/// \param dr       The measured value of the restraint coordinate among its participating atoms
template <typename T>
Vec3<T> restraintDelta(const Vec2<T> init_k, const Vec2<T> final_k, const Vec4<T> init_r,
                       const Vec4<T> final_r, const Vec2<T> mixwt, T dr);

/// \brief Compute the mixture of end-point values that will determine the actual strength and
///        displacement settings of a flat-bottom bimodal harmonic restraint.  The flag about a
///        RestraintApparatus having time-dependent restraints is mostly for convenience, a way to
///        tell whether there is any time-dependent restraint in the collection at all.  Initial
///        and final settings of the steps for each restraint encode whether there is actual time
///        dependence in the result.
///
/// \param step_number  The current step number of the simulation (may include energy minimization
///                     step counts)
/// \param init_step    The initial step at which the restraint engages
/// \param final_step   The final step at which the restraint becomes mature
template <typename T>
Vec2<T> computeRestraintMixture(int step_number, int init_step, int final_step);

} // namespace restraints
} // namespace stormm

#include "restraint_util.tpp"

#endif
