// -*-c++-*-
#ifndef STORMM_TEMP_DISTRIBUTIONS_H
#define STORMM_TEMP_DISTRIBUTIONS_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/symbol_values.h"
#include "Potential/energy_enumerators.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/trajectory_enumerators.h"

namespace stormm {
namespace sampling {

using constants::PrecisionModel;
using constants::ExceptionResponse;
using symbols::boltzmann_constant;
using topology::AtomGraph;

/// \param nw       Number of water molecules
/// \param np       Number of protein atoms
/// \param temp     Temperature at a given instance (float or double)
/// \param f_ener   Energy loss due to constraints
double calc_mu(int n_w, int n_p, double temp, double f_ener);

/// \param m12      The average energy between two given models
/// \param s12      Standard deviation of energies between two given models
/// \param cc       
/// \param u        
double myeval(double m12, double s12, double cc, double u);

/// \param m12      The average energy between two given models
/// \param s12      Standard deviation of energies between two given models
/// \param cc       
double myintegral(double m12, double s12, double cc);

/// \brief Compute the temperature distributions for a given REMD simulation
///
/// \param p_des   Desired Exchange Probability for a swap
/// \param t_low   Lower bound of temperature
/// \param t_high  Higher bound of temperature
/// \param n_w     Number of water molecules
/// \param n_p     Number of protein atoms
/// \param tol     Tolerance for convergence criterion between p_des and the predicted probability
/// \param p_c     Constraints in protein [0: Fully flexible, 1: Bonds to H only, 2: All bonds]
/// \param w_c     Constraints in water [0: Fully flexible, 1: Flexible angle, 2: Rigid]
/// \param h_ff    Hydrogens in Protein [0: All H, 1: Virtual Hydrogen]
/// \param vs      Virtual sites in protein [0: None, 1: Virtual Hydrogen]
/// \param alg     Simulation Type [0: NPT, 1: NVT]
std::vector<double> vanDerSpoel(double p_des, double t_low, double t_high, int n_w, int n_p,
                                double tol, int pc, int wc, int hff, int vs, int alg,
                                ExceptionResponse policy);

} // namespace sampling
} // namespace stormm

#endif
