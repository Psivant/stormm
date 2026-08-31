// -*-c++-*-
#ifndef MMGBSA_NAMELIST_H
#define MMGBSA_NAMELIST_H

#include "../../../src/copyright.h"

namespace permeability {

constexpr int default_simulation_count = 1;
constexpr double default_discard_time = 10.0;
constexpr double default_high_dielectric = 80.0;
constexpr double default_low_dielectric = 4.8;
  
/// \brief Encapsulate the user inputs for a membrane permeability estimator, including its
///        training protocol.
class PermeabilityControls {
public:

private:

  ExceptionResponse policy;   ///< The course of action to take in the event of bad user input
  bool train_model;           ///< Indicate whethere the program will train a model or implement
                              ///<   one based on new (e.g. test) data
  double delta_dipole_wt;     ///< Weight assigned to the change in dipole moment magnitude
                              ///<   observed between the low and high dielectric simulations
  double delta_radgyr_wt;     ///< Weight assigned to the change in radius of gyration observed
                              ///<   between the low and high dielectric simulations
  double delta_
  int simulation_count;       ///< The number of independent simulations to perform for each
                              ///<   molecule in the training set or test set
  double discard_time;        ///< The proportion of each simulation (out of all steps, including
                              ///<   thermal equilibration) to discard before analyzing the
                              ///<   remaining snapshots
  double high_dielectric;     ///< Dielectric constant of the "high dielectric" medium, defaulting
                              ///<   to 80.0 for water
  double low_dielectric;      ///< Dielectric constant of the "low dielectric" medium, defaulting
                              ///<   to 4.0 for protein interiors or the lipid core of a membrane
};

} // namespace permeability

#endif
