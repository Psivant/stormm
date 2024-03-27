// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace topology {

//-------------------------------------------------------------------------------------------------
template <typename T>
VdwCombiningRule inferCombiningRule(const T* lja, const T* ljb, const int lj_type_count,
                                    const ExceptionResponse policy, const bool seek_prevalent) {

  // Trying to determine the combining rule of just one Lennard-Jones parameter is nonsensical.
  // All of the rules could be valid.
  if (lj_type_count <= 1) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("No combining rule can be inferred for " + std::to_string(lj_type_count) +
            " atom types.", "inferCombiningRule");
    case ExceptionResponse::WARN:
      rtWarn("No combining rule can be inferred for " + std::to_string(lj_type_count) +
             " atom types.  A " + getEnumerationName(VdwCombiningRule::LORENTZ_BERTHELOT) +
             " rule will be assumed.", "inferCombiningRule");
      return VdwCombiningRule::LORENTZ_BERTHELOT;
    case ExceptionResponse::SILENT:
      return VdwCombiningRule::LORENTZ_BERTHELOT;
    }
  }
  
  // Check for symmetry in the matrix
  bool is_symmetric = true;
  for (int i = 1; i < lj_type_count; i++) {
    for (int j = 0; j < i; j++) {
      const size_t ij_idx = (i * lj_type_count) + j;
      const size_t ji_idx = (j * lj_type_count) + i;
      is_symmetric = (is_symmetric && fabs(lja[ij_idx] - lja[ji_idx]) < constants::tiny &&
                      fabs(ljb[ij_idx] - ljb[ji_idx]) < constants::tiny);
    }
  }
  if (is_symmetric == false) {
    rtErr("Parameter matrices must be symmetric in order to infer a combining rule",
          "inferCombiningRule");
  }

  // Compute sigma and epsilon parameters based on self interactions.
  std::vector<double> sigma(lj_type_count), epsilon(lj_type_count);
  for (int i = 0; i < lj_type_count; i++) {
    const size_t ii_idx = i * (lj_type_count + 1);
    if (ljb[ii_idx] < constants::small) {
      sigma[i] = 0.0;
      epsilon[i] = 0.0;
    }
    else {
      sigma[i] = sqrt(cbrt(lja[ii_idx] / ljb[ii_idx]));
      double sig_six = sigma[i] * sigma[i] * sigma[i];
      sig_six *= sig_six;
      epsilon[i] = 0.25 * ljb[ii_idx] / sig_six;
    }
  }
  
  // Test for the Lorentz-Berthelot and geometric combining rules
  const std::vector<VdwCombiningRule> trials = { VdwCombiningRule::LORENTZ_BERTHELOT,
                                                 VdwCombiningRule::GEOMETRIC };
  std::vector<int> trial_votes(trials.size(), 0);
  for (size_t i = 0; i < trials.size(); i++) {
    bool test_passes = true;
    int ti = 1;
    while ((test_passes || seek_prevalent) && ti < lj_type_count) {
      int tj = 0;
      while ((test_passes || seek_prevalent) && tj < ti) {
        double pred_eps;
        if (epsilon[ti] < constants::small || epsilon[tj] < constants::small) {
          pred_eps = 0.0;
        }
        else {
          pred_eps = sqrt(epsilon[ti] * epsilon[tj]);
        }
        double pred_sig;
        switch (trials[i]) {
        case VdwCombiningRule::LORENTZ_BERTHELOT:
          pred_sig = 0.5 * (sigma[ti] + sigma[tj]);
          break;
        case VdwCombiningRule::GEOMETRIC:
          if (sigma[ti] < constants::small || sigma[tj] < constants::small) {
            pred_sig = 0.0;
          }
          else {
            pred_sig = sqrt(sigma[ti] * sigma[tj]);
          }
          break;
        case VdwCombiningRule::NBFIX:
          break;
        }
        double pred_sig_six = pred_sig * pred_sig * pred_sig;
        pred_sig_six *= pred_sig_six;
        const double pred_ljb = 4.0 * pred_eps * pred_sig_six;
        const double pred_lja = pred_ljb * pred_sig_six;
        const size_t tij_idx = (tj * lj_type_count) + ti;
        if (seek_prevalent) {
          if (pred_ljb < constants::small) {
            trial_votes[i] += (fabs(pred_lja - lja[tij_idx]) < 1.0e-5 &&
                               fabs(pred_ljb - ljb[tij_idx]) < 1.0e-5);
          }
          else {
            trial_votes[i] += (fabs(pred_lja - lja[tij_idx]) / lja[tij_idx] < 1.0e-5 &&
                               fabs(pred_ljb - ljb[tij_idx]) / ljb[tij_idx] < 1.0e-5);
          }
        }
        else {
          if (pred_ljb < constants::small) {
            test_passes = (test_passes && fabs(pred_lja - lja[tij_idx]) < 1.0e-5 &&
                           fabs(pred_ljb - ljb[tij_idx]) < 1.0e-5);
          }
          else {
            test_passes = (test_passes &&
                           fabs(pred_lja - lja[tij_idx]) / lja[tij_idx] < 1.0e-5 &&
                           fabs(pred_ljb - ljb[tij_idx]) / ljb[tij_idx] < 1.0e-5);
          }
        }
        tj++;
      }
      ti++;
    }
    if (seek_prevalent == false && test_passes) {
      return trials[i];
    }
  }

  // Return the most prevalent combining rule, if there is one of the typical combining rules
  // which can describe at least one off-diagonal interaction and succeeds in describing more
  // such interactions than any other.  If there is a tie, the current list of accepted rules
  // will favor a "Lorentz-Berthelot" characterization.
  if (seek_prevalent) {
    int best_trial_votes = 0;
    int best_trial = -1;
    const int n_trials = trials.size();
    for (int i = 0; i < n_trials; i++) {
      if (trial_votes[i] > best_trial_votes) {
        best_trial_votes = trial_votes[i];
        best_trial = i;
      }
    }
    if (best_trial >= 0) {
      return trials[best_trial];
    }
  }

  // Both known combining rules failed to describe the matrix, and neither made a successful
  // characterization of an interaction between two distinct atom types.
  return VdwCombiningRule::NBFIX;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
VdwCombiningRule inferCombiningRule(const std::vector<T> &lja, const std::vector<T> &ljb,
                                    const ExceptionResponse policy, const bool seek_prevalent) {
  const int n_lj_types = inferLennardJonesTypeCount(lja.size(), ljb.size(), "inferCombiningRule");
  if (n_lj_types * n_lj_types != static_cast<int>(lja.size())) {
    rtErr("Combining rules can only be inferred for square parameter matrices.",
          "inferCombiningRule");
  }
  return inferCombiningRule(lja.data(), ljb.data(), n_lj_types, policy, seek_prevalent);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
VdwCombiningRule inferCombiningRule(const Hybrid<T> &lja, const Hybrid<T> &ljb,
                                    const ExceptionResponse policy, const bool seek_prevalent) {
  const int n_lj_types = inferLennardJonesTypeCount(lja.size(), ljb.size(), "inferCombiningRule");
  return inferCombiningRule(lja.data(), ljb.data(), n_lj_types, policy, seek_prevalent);
}

} // namespace topology
} // namespace stormm
