#include <string>
#include <vector>
#include "copyright.h"
#ifdef STORMM_USE_HPC
#  include "Accelerator/gpu_details.h"
#  include "Accelerator/hpc_config.h"
#endif
#include "Constants/behavior.h"
#include "Constants/symbol_values.h"
#include "DataTypes/stormm_vector_types.h"
#include "FileManagement/file_listing.h"
#include "Math/matrix_ops.h"
#include "Math/vector_ops.h"
#include "MolecularMechanics/dynamics.h"
#include "MolecularMechanics/kinetic.h"
#include "MolecularMechanics/minimization.h"
#include "MolecularMechanics/mm_evaluation.h"
#include "Namelists/nml_minimize.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "Parsing/parsing_enumerators.h"
#include "Potential/energy_enumerators.h"
#include "Potential/scorecard.h"
#include "Random/random.h"
#include "Reporting/error_format.h"
#include "Reporting/summary_file.h"
#include "Structure/local_arrangement.h"
#include "Structure/rattle.h"
#include "Structure/virtual_site_handling.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/integration.h"
#include "Trajectory/motion_sweeper.h"
#include "Trajectory/phasespace.h"
#include "Trajectory/thermostat.h"
#include "Trajectory/trim.h"
#include "UnitTesting/approx.h"
#include "UnitTesting/stopwatch.h"
#include "UnitTesting/test_environment.h"
#include "UnitTesting/test_system_manager.h"
#include "UnitTesting/unit_test.h"

#ifdef STORMM_USE_HPC
using namespace stormm::card;
#endif
using namespace stormm::constants;
using namespace stormm::data_types;
using namespace stormm::diskutil;
using namespace stormm::display;
using namespace stormm::energy;
using namespace stormm::errors;
using namespace stormm::mm;
using namespace stormm::namelist;
using namespace stormm::parse;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::structure;
using namespace stormm::symbols;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Create a pair of point masses tether to one another by a flexible spring or rigid bond.  This
// object can then be simulated as a free object to test energy conservation or under various
// thermostats.
//-------------------------------------------------------------------------------------------------
class Dumbell {
public:

  // The constructor takes parameters for the weights at either end as well as the connection.
  Dumbell(double wa_in = 1.0, double wb_in = 1.0, double keq_in = 1.0, double leq_in = 1.0);

  // Get the mass of the first particle
  double getFirstParticleMass() const;

  // Get the mass of the second particle
  double getSecondParticleMass() const;

  // Get the spring constant connecting the two particles
  double getSpringConstant() const;
  
  // Get the equilibrium length of the connection between the two particles
  double getConnectionLength() const;
  
  // Set the mass of the first particle
  void setFirstParticleMass(double wa_in);

  // Set the mass of the second particle
  void setSecondParticleMass(double wb_in);

  // Set the spring constant connecting the two particles
  void setSpringConstant(double keq_in);

  // Set the equilibrium length of the bond between the two particles
  void setConnectionLength(double leq_in);

  // Compute the force acting on the two particles
  //
  // Arguments:
  //   r:  The positions of each particle
  //   f:  The forces acting on each particle, computed and returned
  void computeForce(const std::vector<double3> &r, std::vector<double3> *f) const;

  // Compute the kinetic and potential energy of the system, in kcal/mol.  The result is returned
  // as a tuple with the kinetic energy in its "x" member and the potential energy in its "y"
  // member.
  //
  // Arguments:
  //   r:  The positions of each particle
  //   v:  The velocities of each particle
  double2 computeEnergy(const std::vector<double3> &r, const std::vector<double3> &v) const;

  // Compute the net system momentum of the dumbell.
  //
  // Arguments:
  //   v:  The velocities of each particle
  double3 computeNetMomentum(const std::vector<double3> &v) const;

  // Update the particle velocities based on computed forces.  This routine works in the context of
  // a velocity-Verlet integration step, updating velocities according to half the time step.  Call
  // it twice to update velocities over the full course of the time step.
  //
  // Arguments:
  //   f:   Forces on each particle
  //   v:   Velocities of each particle, updated and returned
  //   dt:  The time step 
  void updateVelocities(const std::vector<double3> &f, std::vector<double3> *v, double dt) const;

  // Update the particle positions based on their current velocities.  This routine works in the
  // context of a velocity-Verlet integration step, updating positions according to the full time
  // step based on the velocities computed after half the impulse due to the current forces has
  // been contributed.
  //
  // Arguments:
  //   v:   Velocities of each particle
  //   r:   Positions of each particle, updated and returned
  //   dt:  The time step
  void updatePositions(const std::vector<double3> &v, std::vector<double3> *r, double dt) const;
  
private:
  double wa;   // The mass of the first weight, in atomic units
  double wb;   // The mass of the second weight, in atomic units
  double keq;  // The stiffness of the connecting spring (if a rigid bond approximation is used,
               //   this parameter will have no effect), in kcal/mol-A^2
  double leq;  // The equilibrium length of the connecting spring, in Angstroms
};

//-------------------------------------------------------------------------------------------------
Dumbell::Dumbell(const double wa_in, const double wb_in, const double keq_in,
                 const double leq_in) :
  wa{wa_in}, wb{wb_in}, keq{keq_in}, leq{leq_in}
{}

//-------------------------------------------------------------------------------------------------
double Dumbell::getFirstParticleMass() const {
  return wa;
}

//-------------------------------------------------------------------------------------------------
double Dumbell::getSecondParticleMass() const {
  return wb;
}

//-------------------------------------------------------------------------------------------------
double Dumbell::getSpringConstant() const {
  return keq;
}

//-------------------------------------------------------------------------------------------------
double Dumbell::getConnectionLength() const {
  return leq;
}

//-------------------------------------------------------------------------------------------------
void Dumbell::setFirstParticleMass(const double wa_in) {
  wa = wa_in;
}

//-------------------------------------------------------------------------------------------------
void Dumbell::setSecondParticleMass(const double wb_in) {
  wb = wb_in;
}

//-------------------------------------------------------------------------------------------------
void Dumbell::setSpringConstant(const double keq_in) {
  keq = keq_in;
}

//-------------------------------------------------------------------------------------------------
void Dumbell::setConnectionLength(const double leq_in) {
  leq = leq_in;
}

//-------------------------------------------------------------------------------------------------
void Dumbell::computeForce(const std::vector<double3> &r, std::vector<double3> *f) const {
  const double3 dr = { r[1].x - r[0].x, r[1].y - r[0].y, r[1].z - r[0].z };
  const double lnew = sqrt((dr.x * dr.x) + (dr.y * dr.y) + (dr.z * dr.z));
  const double dl = (lnew - leq);
  const double fmag = (lnew < tiny) ? 2.0 * keq * dl / tiny : 2.0 * keq * dl / lnew;
  const double fx = fmag * dr.x;
  const double fy = fmag * dr.y;
  const double fz = fmag * dr.z;
  double3* f_ptr = f->data();
  f_ptr[0] = { fx, fy, fz };
  f_ptr[1] = { -fx, -fy, -fz };
}

//-------------------------------------------------------------------------------------------------
double2 Dumbell::computeEnergy(const std::vector<double3> &r,
                               const std::vector<double3> &v) const {
  double2 result;

  // Compute the kinetic energy
  const double vasq = (v[0].x * v[0].x) + (v[0].y * v[0].y) + (v[0].z * v[0].z);
  const double vbsq = (v[1].x * v[1].x) + (v[1].y * v[1].y) + (v[1].z * v[1].z);
  result.x = 0.5 * gafs_to_kcal * ((wa * vasq) + (wb * vbsq));
  
  // Compute the potential energy
  const double3 dr = { r[1].x - r[0].x, r[1].y - r[0].y, r[1].z - r[0].z };
  const double lnew = sqrt((dr.x * dr.x) + (dr.y * dr.y) + (dr.z * dr.z));
  const double dl = (lnew - leq);
  result.y = keq * dl * dl;
  return result;
}

//-------------------------------------------------------------------------------------------------
double3 Dumbell::computeNetMomentum(const std::vector<double3> &v) const {
  double3 result = { (wa * v[0].x) + (wb * v[1].x),
                     (wa * v[0].y) + (wb * v[1].y),
                     (wa * v[0].z) + (wb * v[1].z) };
  return result;
}

//-------------------------------------------------------------------------------------------------
void Dumbell::updateVelocities(const std::vector<double3> &f, std::vector<double3> *v,
                               const double dt) const {
  double3* v_ptr = v->data();
  const double hmadt = 0.5 * dt * kcal_to_gafs / wa;
  v_ptr[0].x += hmadt * f[0].x;
  v_ptr[0].y += hmadt * f[0].y;
  v_ptr[0].z += hmadt * f[0].z;
  const double hmbdt = 0.5 * dt * kcal_to_gafs / wb;
  v_ptr[1].x += hmbdt * f[1].x;
  v_ptr[1].y += hmbdt * f[1].y;
  v_ptr[1].z += hmbdt * f[1].z;
}

//-------------------------------------------------------------------------------------------------
void Dumbell::updatePositions(const std::vector<double3> &v, std::vector<double3> *r,
                              const double dt) const {

  // Position updates are performed once per step.
  double3* r_ptr= r->data();
  r_ptr[0].x += dt * v[0].x;
  r_ptr[0].y += dt * v[0].y;
  r_ptr[0].z += dt * v[0].z;
  r_ptr[1].x += dt * v[1].x;
  r_ptr[1].y += dt * v[1].y;
  r_ptr[1].z += dt * v[1].z;
}

//-------------------------------------------------------------------------------------------------
// Set up a series of test systems to assess momentum removal.  Begin with one that has a net
// translational momentum and a displacement from the origin, but no rotation.  Next, try one with
// rotational momentum along the Cartesian z-axis only.
//-------------------------------------------------------------------------------------------------
void runMomentumTests() {

  // Set up a system with only translational momentum.
  std::vector<double> xcrd(7, 0.0), ycrd(7, 0.0), zcrd(7, 0.0);
  std::vector<double> xvel(7, 0.0), yvel(7, 0.0), zvel(7, 0.0);
  std::vector<double> masses(7, 1.0);
  xcrd[1] =  1.0;
  xcrd[2] = -1.0;
  ycrd[3] =  1.0;
  ycrd[4] = -1.0;
  zcrd[5] =  1.0;
  zcrd[6] = -1.0;
  const std::vector<double> first_x_ans = xcrd;
  const std::vector<double> first_y_ans = ycrd;
  const std::vector<double> first_z_ans = zcrd;
  for (int i = 0; i < 7; i++) {
    xcrd[i] += 1.1;
    ycrd[i] -= 0.7;
    zcrd[i] += 3.6;
  }
  xvel[0] = -6.1;
  yvel[0] =  2.7;
  zvel[0] = -1.1;
  const double mean_xv = mean(xvel);
  const double mean_yv = mean(yvel);
  const double mean_zv = mean(zvel);
  std::vector<double> first_xv_ans(7), first_yv_ans(7), first_zv_ans(7);
  for (int i = 0; i < 7; i++) {
    first_xv_ans[i] = xvel[i] - mean_xv;
    first_yv_ans[i] = yvel[i] - mean_yv;
    first_zv_ans[i] = zvel[i] - mean_zv;
  }
  removeMomentum<double, double, double>(xcrd.data(), ycrd.data(), zcrd.data(), nullptr, nullptr,
                                         nullptr, xvel.data(), yvel.data(), zvel.data(), nullptr,
                                         nullptr, nullptr, masses.data(), UnitCellType::NONE, 7);
  check(xcrd, RelationalOperator::EQUAL, first_x_ans, "Re-centering along the X axis fails on a "
        "system with only net translational momentum.");
  check(ycrd, RelationalOperator::EQUAL, first_y_ans, "Re-centering along the Y axis fails on a "
        "system with only net translational momentum.");
  check(zcrd, RelationalOperator::EQUAL, first_z_ans, "Re-centering along the Z axis fails on a "
        "system with only net translational momentum.");
  check(xvel, RelationalOperator::EQUAL, first_xv_ans, "Translational momentum purge along the X "
        "axis fails in the simplest system.");
  check(yvel, RelationalOperator::EQUAL, first_yv_ans, "Translational momentum purge along the Y "
        "axis fails in the simplest system.");
  check(zvel, RelationalOperator::EQUAL, first_zv_ans, "Translational momentum purge along the Z "
        "axis fails in the simplest system.");

  // Add net rotational momentum to the system along the Cartesian Z-axis
  for (int i = 0; i < 7; i++) {
    xvel[i] = 0.0;
    yvel[i] = 0.0;
    zvel[i] = 0.0;
  }
  yvel[1] =  1.0;
  yvel[2] = -1.0;
  removeMomentum<double, double, double>(xcrd.data(), ycrd.data(), zcrd.data(), nullptr, nullptr,
                                         nullptr, xvel.data(), yvel.data(), zvel.data(), nullptr,
                                         nullptr, nullptr, masses.data(), UnitCellType::NONE, 7);
  
  // In this example, the distal points at +1 and -1 along the X-axis rotate counter-clockwise
  // as would been seen looking down the Z-axis.  The solution is to make the entire object spin
  // against that rotation at half the speed (there are four points of equal mass which affect the
  // angular moment of intertia along the Z-axis).  The distal points along the X-axis will
  // continue to move counter-clockwise but only at half their original rate.  Meanwhile, the
  // distal points along the Y-axis will rotate clockwise to counterbalance the angular momentum of
  // te system.
  const std::vector<double> second_xv_ans = {  0.0,  0.0,  0.0,  0.5, -0.5,  0.0,  0.0 };
  const std::vector<double> second_yv_ans = {  0.0,  0.5, -0.5,  0.0,  0.0,  0.0,  0.0 };
  const std::vector<double> second_zv_ans = {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 };
  
  check(xvel, RelationalOperator::EQUAL, second_xv_ans, "Removal of rotational momentum from a "
        "simple system fails in a symmetric case.");
  check(yvel, RelationalOperator::EQUAL, second_yv_ans, "Removal of rotational momentum from a "
        "simple system fails in a symmetric case.");
  check(zvel, RelationalOperator::EQUAL, second_zv_ans, "Removal of rotational momentum from a "
        "simple system fails in a symmetric case.");

  // Change the layout such that the particles are now spread differently along each axis.
  xcrd[0] =  0.25;
  xcrd[1] =  0.75;
  ycrd[1] = -0.85;
  xcrd[2] =  0.71;
  xcrd[3] = -0.66;
  for (int i = 0; i < 7; i++) {
    xvel[i] = 0.0;
    yvel[i] = 0.0;
    zvel[i] = 0.0;
  }
  yvel[1] =  1.0;
  yvel[2] = -1.0;
  removeMomentum<double, double, double>(xcrd.data(), ycrd.data(), zcrd.data(), nullptr, nullptr,
                                         nullptr, xvel.data(), yvel.data(), zvel.data(), nullptr,
                                         nullptr, nullptr, masses.data(), UnitCellType::NONE, 7);

  // The angular momentum of the asymmetric arrangement will be messier, but rotation is still
  // confined to the Cartesian Z-axis.
  const std::vector<double> third_xv_ans = {  0.00120630, -0.00723779,  0.00120630,  0.01114052,
                                             -0.00872792,  0.00120630,  0.00120630 };
  const std::vector<double> third_yv_ans = { -0.00099342,  0.99403947, -1.00556316,  0.00804672,
                                              0.00149013,  0.00149013,  0.00149013 };
  const std::vector<double> third_zv_ans(7, 0.0);
  check(xvel, RelationalOperator::EQUAL, third_xv_ans, "Removal of rotational momentum from an "
        "asymmetric system fails with deviations along the Cartesian X-axis.");
  check(yvel, RelationalOperator::EQUAL, third_yv_ans, "Removal of rotational momentum from an "
        "asymmetric system fails with deviations along the Cartesian Y-axis.");
  check(zvel, RelationalOperator::EQUAL, third_zv_ans, "Removal of rotational momentum from an "
        "asymmetric system fails with deviations along the Cartesian Z-axis.");

  // Check other axes
  for (int i = 0; i < 7; i++) {
    xcrd[i] = 0.0;
    ycrd[i] = 0.0;
    zcrd[i] = 0.0;
    xvel[i] = 0.0;
    yvel[i] = 0.0;
    zvel[i] = 0.0;
  }
  xcrd[1] =  1.0;
  xcrd[2] = -1.0;
  ycrd[3] =  0.15;
  ycrd[4] = -1.4;
  zcrd[5] =  0.8;
  zcrd[6] = -2.27;
  zvel[3] = -1.0;
  zvel[4] =  1.0;
  removeMomentum<double, double, double>(xcrd.data(), ycrd.data(), zcrd.data(), nullptr, nullptr,
                                         nullptr, xvel.data(), yvel.data(), zvel.data(), nullptr,
                                         nullptr, nullptr, masses.data(), UnitCellType::NONE, 7);
  const std::vector<double> fourth_xv_ans(7, 0.0);
  const std::vector<double> fourth_yv_ans = { -0.04493693, -0.04493693, -0.04493693, -0.04493693,
                                              -0.04493693, -0.21612523,  0.44080987 };
  const std::vector<double> fourth_zv_ans = {  0.03821167,  0.03821167,  0.03821167, -0.92969052,
                                               0.73863215,  0.03821167,  0.03821167 };
  check(xvel, RelationalOperator::EQUAL, fourth_xv_ans, "Removal of rotational momentum from an "
        "asymmetric system fails with deviations along the Cartesian X-axis.");
  check(yvel, RelationalOperator::EQUAL, fourth_yv_ans, "Removal of rotational momentum from an "
        "asymmetric system fails with deviations along the Cartesian Y-axis.");
  check(zvel, RelationalOperator::EQUAL, fourth_zv_ans, "Removal of rotational momentum from an "
        "asymmetric system fails with deviations along the Cartesian Z-axis.");

  // Make perturbations along all three axes simultaneously
  Xoroshiro128pGenerator xrs;
  for (int i = 0; i < 7; i++) {
    xcrd[i] = 10.0 * xrs.uniformRandomNumber();
    ycrd[i] = 10.0 * xrs.uniformRandomNumber();
    zcrd[i] = 10.0 * xrs.uniformRandomNumber();
    xvel[i] = 0.1 * (xrs.uniformRandomNumber() - 0.5);
    yvel[i] = 0.1 * (xrs.uniformRandomNumber() - 0.5);
    zvel[i] = 0.1 * (xrs.uniformRandomNumber() - 0.5);
  }
  removeMomentum<double, double, double>(xcrd.data(), ycrd.data(), zcrd.data(), nullptr, nullptr,
                                         nullptr, xvel.data(), yvel.data(), zvel.data(), nullptr,
                                         nullptr, nullptr, masses.data(), UnitCellType::NONE, 7);
  const std::vector<double> fifth_xv_ans = { -0.03105789,  0.01654014,  0.00939234,  0.01446993,
                                             -0.01064392, -0.00484942,  0.00614883 };
  const std::vector<double> fifth_yv_ans = { -0.02735157, -0.02311428, -0.02610229,  0.00314351,
                                              0.02743904,  0.04356921,  0.00241638 };
  const std::vector<double> fifth_zv_ans = {  0.02940335,  0.02417684,  0.01718409,  0.01611805,
                                             -0.02733724, -0.02234826, -0.03719683 };
  check(xvel, RelationalOperator::EQUAL, fifth_xv_ans, "Removal of rotational momentum from a "
        "randomized system fails with deviations along the Cartesian X-axis.");
  check(yvel, RelationalOperator::EQUAL, fifth_yv_ans, "Removal of rotational momentum from a "
        "randomized system fails with deviations along the Cartesian Y-axis.");
  check(zvel, RelationalOperator::EQUAL, fifth_zv_ans, "Removal of rotational momentum from a "
        "randomized system fails with deviations along the Cartesian Z-axis.");
}

//-------------------------------------------------------------------------------------------------
// Perform energy-conserving dynamics on a system of two particles connected by a flexible spring.
// Return the result in terms of the mean and standard deviation of the kinetic and potential
// energy, as output by the Dumbell class's computeEnergy() member function.
//
// Arguments:
//   bar:            The dumbell to simulate
//   dt:             The time step to use, in femtoseconds
//   nstep:          The length of simulation to perform
//   blk_ave_step:   Number of steps over which to take block averages for the resulting energies
//   init_r:         Initial positions of the particles, in Angstroms
//   init_v:         Initial velocities of the particles, in Angstroms per femtosecond
//   total_nrg_var:  Variance of the total energy over the entire simulation (computed and
//                   returned, if the pointer is valid)
//-------------------------------------------------------------------------------------------------
std::vector<double4> isoenergetic(const Dumbell &bar, const double dt, const int nstep,
                                  const int blk_ave_step, const std::vector<double3> &init_r,
                                  const std::vector<double3> &init_v,
                                  double *total_nrg_var = nullptr) {
  std::vector<double3> r(2), r_ref(2), v(2), v_ref(2), f(2);
  const bool calc_e_var = (total_nrg_var != nullptr);

  // Initialize the positions and velocities
  for (int i = 0; i < 2; i++) {
    r[i] = init_r[i];
    v[i] = init_v[i];
  }
  
  // Run the simulation
  std::vector<double4> result;
  result.reserve((nstep / blk_ave_step) - 1);
  int blk_idx = 0;
  double2 blk_ave, blk_ave_sq;
  const double bstep = blk_ave_step;
  const double inv_bstep = 1.0 / static_cast<double>(blk_ave_step);
  double total_e = 0.0;
  double total_e_sq = 0.0;
  for (int nt = 0; nt < nstep; nt++) {

    // Compute the force on the two particles
    bar.computeForce(r, &f);

    // Velocity-Verlet 1v: Update the velocities by half a step
    bar.updateVelocities(f, &v, dt);

    // Compute the energy
    const double2 u = bar.computeEnergy(r, v);
    if (nt % blk_ave_step == 0) {
      if (nt > blk_ave_step) {
        const double ke_ave = blk_ave.x * inv_bstep;
        const double pe_ave = blk_ave.y * inv_bstep;
        const double ke_std = sqrt((bstep * blk_ave_sq.x) - (blk_ave.x * blk_ave.x)) /
                              sqrt(bstep * (bstep - 1.0));
        const double pe_std = sqrt((bstep * blk_ave_sq.y) - (blk_ave.y * blk_ave.y)) /
                              sqrt(bstep * (bstep - 1.0));
        result.push_back({ ke_ave, pe_ave, ke_std, pe_std });
      }
      blk_ave = { 0.0, 0.0 };
      blk_ave_sq = { 0.0, 0.0 };
    }
    blk_ave.x += u.x;
    blk_ave.y += u.y;
    blk_ave_sq.x += u.x * u.x;
    blk_ave_sq.y += u.y * u.y;
    if (calc_e_var) {
      total_e += u.x + u.y;
      total_e_sq += (u.x + u.y) * (u.x + u.y);
    }
    
    // Velocity-Verlet 2v: Update the velocities by half a step
    bar.updateVelocities(f, &v, dt);

    // Velocity-Verlet 1c: Update the positions
    bar.updatePositions(v, &r, dt);
  }
  if (calc_e_var) {
    const double dnstep = nstep;
    *total_nrg_var = sqrt((dnstep * total_e_sq) - (total_e * total_e)) /
                     sqrt(dnstep * (dnstep - 1.0));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Carry out the MD cycle in vacuum phase.  The precise origin (orientation) of the PhaseSpace
// abstract can be toggled to work from one point in the time cycle and construct the conditions
// for the next.
//
// Arguments:
//   
//-------------------------------------------------------------------------------------------------
void mockDynamicsStep(PhaseSpaceWriter *psw, ScoreCard *sc, const ThermostatReader<double> &tstr,
                      const ValenceKit<double> &vk, const NonbondedKit<double> &nbk,
                      const StaticExclusionMaskReader &ser, const VirtualSiteKit<double> &vsk,
                      const ConstraintKit<double> &cnk, const double dt, const double cnst_tol,
                      const RattleMethod rattle_strat, const ChemicalDetailsKit &cdk) {

  // Compute the forces on all atoms due to all vacuum-phase MM terms
  evalNonbValeMM(*psw, sc, vk, nbk, ser, EvaluateForce::YES, 0, 0.0, 0.0);
  sc->commit(0);

  // Transfer virtual site forces
  transmitVirtualSiteForces(*psw, vsk);
  
  // Advance particle velocities.
  velocityVerletVelocityUpdate(psw, cdk, tstr);

  // Apply velocity constraints
  if (cnst_tol > 0.0) {
    rattleVelocities<double>(psw, cnk, dt, cnst_tol, 20, rattle_strat);
  }

  // Compute the kinetic energy.
  evalKineticEnergy(*psw, sc, cdk, 0);

  // Advance particle velocities once more, then advance particle positions.
  velocityVerletCoordinateUpdate(psw, cdk, tstr);

  // Apply positional constraints
  if (cnst_tol > 0.0) {
    rattlePositions<double>(psw, cnk, dt, cnst_tol, 20, rattle_strat);
  }

  // Place virtual sites in the new coordinates
  placeVirtualSites<double, double>(psw->xalt, psw->yalt, psw->zalt, psw->umat_alt, psw->invu_alt,
                                    psw->unit_cell, vsk);
  
  // Re-initialize force accumulators in preparation for the next step
  for (int i = 0; i < psw->natom; i++) {
    psw->fxalt[i] = 0.0;
    psw->fyalt[i] = 0.0;
    psw->fzalt[i] = 0.0;
  }
}

//-------------------------------------------------------------------------------------------------
// Carry out a simple molecular dynamics simulation of a small molecule, with no thermostating.
// Return the mean temperature from the run.
//
// Arguments:
//   ps:            Coordinates, velocities, and forces of all particles in the system
//   ag:            Topology for the system, containing all parameters
//   nstep:         The number of steps to take
//   dt:            The time step, in femtoseconds
//   rattle_tol:    The RATTLE tolerance (set to zero to suppress bond constraint applications)
//   rattle_strat:  The manner by which to iterate on hub-and-spoke constraint groups
//   do_tests:      Indicate whether tets are possible
//   nblock:        Number of time steps to group into each block for block averaging
//   equipart_tol:  Criterion for valid equipartition
//   nrg_cons_tol:  Criterion for successful energy conservation
//-------------------------------------------------------------------------------------------------
double mdTrial(PhaseSpace *ps, const AtomGraph &ag, const int nstep, const double dt,
               const double rattle_tol, const RattleMethod rattle_strat,
               const TestPriority do_tests, const int nblock = 1000,
               const double equipart_tol = 3.0e-2, const double nrg_cons_tol = 4.0e-3) {
  double result = 0.0;
  
  // Create a thermostat to control the integration
  Thermostat tst(ag, ThermostatKind::NONE);
  tst.setTimeStep(dt);
  ThermostatWriter<double> tstw = tst.dpData();
  const ThermostatReader<double> tstr(tstw);

  // Create an energy tracker
  ScoreCard sc(1, nstep / nblock, 32);

  // Obtain critical abstracts of the topology
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
  const StaticExclusionMask se(ag);
  const StaticExclusionMaskReader ser = se.data();
  const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
  const ConstraintKit<double> cnk = ag.getDoublePrecisionConstraintKit();
  
  // Assume that the valid coordinates are present at the WHITE point in the coordinate cycle.
  // While each abstract sees positions, velocities, and forces at both points in the time cycle,
  // make two abstracts so that toggling between them ("tick-tock") always phrases the force
  // calculation in terms of the "WHITE" point in the coordinate cycle.
  PhaseSpaceWriter psw = ps->data(CoordinateCycle::WHITE);
  PhaseSpaceWriter psw_alt = ps->data(CoordinateCycle::BLACK);

  // Adjust the virtual sites to the exact coordinates (do not trust the original file, which may
  // not have a very high precision format).
  placeVirtualSites(psw, vsk);

  // Iterate through all requested dynamics steps
  double pe = 0.0;
  double ke = 0.0;
  double te = 0.0;
  double pe_sq = 0.0;
  double ke_sq = 0.0;
  double te_sq = 0.0;
  const int actual_nstep = roundUp(nstep, nblock);
  const int npts = actual_nstep / nblock;
  const double dnblock = nblock;
  std::vector<double> mean_pe(npts), mean_ke(npts), mean_te(npts);
  std::vector<double> stdev_pe(npts), stdev_ke(npts), stdev_te(npts);
  for (int i = 0; i < actual_nstep; i++) {

    // Run the dynamics step
    if (i & 0x1) {
      mockDynamicsStep(&psw_alt, &sc, tstr, vk, nbk, ser, vsk, cnk, dt, rattle_tol, rattle_strat,
                       cdk);
    }
    else {
      mockDynamicsStep(&psw, &sc, tstr, vk, nbk, ser, vsk, cnk, dt, rattle_tol, rattle_strat, cdk);
    }

    // Udpate the coordiante object's cycle position
    ps->updateCyclePosition();
    
    const double dpe = sc.reportPotentialEnergy(0);
    const double dke = sc.reportInstantaneousStates(StateVariable::KINETIC, 0);

    // Compute the temperature and accumulate the result
    computeTemperature(&sc, ag, true, 0);
    result += sc.reportInstantaneousStates(StateVariable::TEMPERATURE_ALL, 0);

    // Accumulate potential and kinetic energies.  Track the variance of each.
    pe += dpe;
    ke += dke;
    te += dpe + dke;
    pe_sq += dpe * dpe;
    ke_sq += dke * dke;
    te_sq += (dpe + dke) * (dpe + dke);
    if ((i + 1) % nblock == 0) {
      const size_t pt_idx = i / nblock;
      mean_pe[pt_idx] = pe / dnblock;
      mean_ke[pt_idx] = ke / dnblock;
      mean_te[pt_idx] = te / dnblock;
      stdev_pe[pt_idx] = sqrt((dnblock * pe_sq) - (pe * pe)) /
                         sqrt(dnblock * (dnblock - 1.0));
      stdev_ke[pt_idx] = sqrt((dnblock * ke_sq) - (ke * ke)) /
                         sqrt(dnblock * (dnblock - 1.0));
      stdev_te[pt_idx] = sqrt((dnblock * te_sq) - (te * te)) /
                         sqrt(dnblock * (dnblock - 1.0));
      pe = 0.0;
      ke = 0.0;
      te = 0.0;
      pe_sq = 0.0;
      ke_sq = 0.0;
      te_sq = 0.0;
    }
  }
  
  // Test for equipartition
  check(mean(stdev_pe), RelationalOperator::EQUAL,
        Approx(mean(stdev_ke)).margin(equipart_tol * dt), "Equipartition is violated for a "
        "simulation of " + getBaseName(ag.getFileName()) + " with a time step of " +
        minimalRealFormat(dt, 1.0e-4) + " fs.", do_tests);

  // Test for energy conservation
  if (npts > 2) {
    check(mean_te[npts - 1], RelationalOperator::EQUAL,
          Approx(mean_te[1]).margin(nrg_cons_tol * dt), "A simulation of " +
          getBaseName(ag.getFileName()) + " with a time step of " + minimalRealFormat(dt, 1.0e-4) +
          "fs does not appear to conserve energy.", do_tests);
  }

  return result / static_cast<double>(actual_nstep);
}

//-------------------------------------------------------------------------------------------------
// Check the removal of momentum from systems in a synthesis, using the MotionSweeper class built
// for calculations on the GPU.
//
// Arguments:
//   tsm:  A collection of test systems that will form the synthesis
//   xrs:  Source of random numbers to provide instabilities that the momentum removal will correct
//   gpu:  Specifications of the GPU that will make use of the MotionSweeper.  The default value of
//         a null GPU will work entirely off of the class's CPU routines.
//-------------------------------------------------------------------------------------------------
void testMotionSweeper(const TestSystemManager &tsm, Xoshiro256ppGenerator *xrs,
                       const GpuDetails &gpu = null_gpu) {
  HybridTargetLevel insp_tier;
  std::string gpu_engaged_msg("");
  if (gpu == null_gpu) {
    insp_tier = HybridTargetLevel::HOST;
  }
#ifdef STORMM_USE_HPC
  else {
    gpu_engaged_msg = "  The GPU kernels were engaged.";
    insp_tier = HybridTargetLevel::DEVICE;
  }
#endif
  const std::vector<int> mover_idx = incrementingSeries(0, tsm.getSystemCount());
  PhaseSpaceSynthesis mover_ps = tsm.exportPhaseSpaceSynthesis(mover_idx);
  AtomGraphSynthesis mover_ag = tsm.exportAtomGraphSynthesis(mover_idx);
  MotionSweeper mos(mover_ps);
  PsSynthesisWriter mover_psw = mover_ps.data(insp_tier);
  PsSynthesisWriter host_mover_psw = mover_ps.data();
  MotionSweepWriter mosw = mos.data(insp_tier);
  const SyAtomUpdateKit<double,
                        double2,
                        double4> mover_auk = mover_ag.getDoublePrecisionAtomUpdateKit(insp_tier);
  std::vector<std::vector<double>> x_prtb, y_prtb, z_prtb, xv_prtb, yv_prtb, zv_prtb;
  x_prtb.resize(host_mover_psw.system_count);
  y_prtb.resize(host_mover_psw.system_count);
  z_prtb.resize(host_mover_psw.system_count);
  xv_prtb.resize(host_mover_psw.system_count);
  yv_prtb.resize(host_mover_psw.system_count);
  zv_prtb.resize(host_mover_psw.system_count);
  for (int i = 0; i < host_mover_psw.system_count; i++) {
    x_prtb[i] = gaussianRand(xrs, host_mover_psw.atom_counts[i], 0.15);
    y_prtb[i] = gaussianRand(xrs, host_mover_psw.atom_counts[i], 0.15);
    z_prtb[i] = gaussianRand(xrs, host_mover_psw.atom_counts[i], 0.15);
    std::vector<double> xyz_i = mover_ps.getInterlacedCoordinates(i);
    for (int j = 0; j < host_mover_psw.atom_counts[i]; j++) {
      x_prtb[i][j] += xyz_i[(3 * j)    ];
      y_prtb[i][j] += xyz_i[(3 * j) + 1];
      z_prtb[i][j] += xyz_i[(3 * j) + 2];
    }
    xv_prtb[i] = gaussianRand(xrs, host_mover_psw.atom_counts[i], 0.1);
    yv_prtb[i] = gaussianRand(xrs, host_mover_psw.atom_counts[i], 0.1);
    zv_prtb[i] = gaussianRand(xrs, host_mover_psw.atom_counts[i], 0.1);
    addScalarToVector<double>(&x_prtb[i],  2.917);
    addScalarToVector<double>(&y_prtb[i], -5.882);
    addScalarToVector<double>(&z_prtb[i],  1.079);
    addScalarToVector<double>(&xv_prtb[i],  0.017);
    addScalarToVector<double>(&yv_prtb[i], -0.009);
    addScalarToVector<double>(&zv_prtb[i], -0.011);
    for (int j = 0; j < host_mover_psw.atom_counts[i]; j++) {
      const int95_t ix = hostDoubleToInt95(x_prtb[i][j] * host_mover_psw.gpos_scale);
      const int95_t iy = hostDoubleToInt95(y_prtb[i][j] * host_mover_psw.gpos_scale);
      const int95_t iz = hostDoubleToInt95(z_prtb[i][j] * host_mover_psw.gpos_scale);
      const int95_t ivx = hostDoubleToInt95(xv_prtb[i][j] * host_mover_psw.vel_scale);
      const int95_t ivy = hostDoubleToInt95(yv_prtb[i][j] * host_mover_psw.vel_scale);
      const int95_t ivz = hostDoubleToInt95(zv_prtb[i][j] * host_mover_psw.vel_scale);
      const size_t atom_idx = j + host_mover_psw.atom_starts[i];
      host_mover_psw.xcrd[atom_idx] = ix.x;
      host_mover_psw.ycrd[atom_idx] = iy.x;
      host_mover_psw.zcrd[atom_idx] = iz.x;
      host_mover_psw.xcrd_ovrf[atom_idx] = ix.y;
      host_mover_psw.ycrd_ovrf[atom_idx] = iy.y;
      host_mover_psw.zcrd_ovrf[atom_idx] = iz.y;
      host_mover_psw.xvel[atom_idx] = ivx.x;
      host_mover_psw.yvel[atom_idx] = ivy.x;
      host_mover_psw.zvel[atom_idx] = ivz.x;
      host_mover_psw.xvel_ovrf[atom_idx] = ivx.y;
      host_mover_psw.yvel_ovrf[atom_idx] = ivy.y;
      host_mover_psw.zvel_ovrf[atom_idx] = ivz.y;
    }
  }

  // Upload the coordinate synthesis in order to make a new copy of the data and prepare to
  // evaluate the GPU kernels.
#ifdef STORMM_USE_HPC
  if (gpu != null_gpu) {
    mover_ag.upload();    
    mover_ps.upload();
    mos.uploadAll();
  }
#endif

  // Compute the centers of mass
  const int nsys = host_mover_psw.system_count;
  accumulateCenterOfMassMotion(&mosw, mover_auk, mover_psw, gpu);
  std::vector<double> synth_comx(nsys), synth_comy(nsys), synth_comz(nsys), tmass(nsys);
  std::vector<double> synth_comx_ans(nsys), synth_comy_ans(nsys), synth_comz_ans(nsys);
  std::vector<double> synth_velx(nsys), synth_vely(nsys), synth_velz(nsys);
  std::vector<double> synth_velx_ans(nsys), synth_vely_ans(nsys), synth_velz_ans(nsys);
  std::vector<std::vector<double>> synth_masses(nsys);
  for (int i = 0; i < nsys; i++) {

    // Center of mass
    const double3 com_i = mos.getCenterOfMass(i, insp_tier);
    synth_comx[i] = com_i.x;
    synth_comy[i] = com_i.y;
    synth_comz[i] = com_i.z;
    synth_comx_ans[i] = 0.0;
    synth_comy_ans[i] = 0.0;
    synth_comz_ans[i] = 0.0;
    const AtomGraph *i_ag = mover_ps.getSystemTopologyPointer(i);
    synth_masses[i] = i_ag->getAtomicMass<double>();
    for (int j = 0; j < host_mover_psw.atom_counts[i]; j++) {
      synth_comx_ans[i] += x_prtb[i][j] * synth_masses[i][j];
      synth_comy_ans[i] += y_prtb[i][j] * synth_masses[i][j];
      synth_comz_ans[i] += z_prtb[i][j] * synth_masses[i][j];
    }
    tmass[i] = sum<double>(synth_masses[i]);
    synth_comx_ans[i] /= tmass[i];
    synth_comy_ans[i] /= tmass[i];
    synth_comz_ans[i] /= tmass[i];

    // Net velocity
    const double3 synth_vi = mos.getNetVelocity(i, insp_tier);
    synth_velx[i] = synth_vi.x;
    synth_vely[i] = synth_vi.y;
    synth_velz[i] = synth_vi.z;
    synth_velx_ans[i] = 0.0;
    synth_vely_ans[i] = 0.0;
    synth_velz_ans[i] = 0.0;
    for (int j = 0; j < host_mover_psw.atom_counts[i]; j++) {
      synth_velx_ans[i] += xv_prtb[i][j] * synth_masses[i][j];
      synth_vely_ans[i] += yv_prtb[i][j] * synth_masses[i][j];
      synth_velz_ans[i] += zv_prtb[i][j] * synth_masses[i][j];
    }
    synth_velx_ans[i] /= tmass[i];
    synth_vely_ans[i] /= tmass[i];
    synth_velz_ans[i] /= tmass[i];
  }
  check(synth_comx, RelationalOperator::EQUAL, synth_comx_ans, "Centers of mass computed using a "
        "MotionSweeper do not agree with basic, independent calculations along the X dimension." +
        gpu_engaged_msg, tsm.getTestingStatus());
  check(synth_comy, RelationalOperator::EQUAL, synth_comy_ans, "Centers of mass computed using a "
        "MotionSweeper do not agree with basic, independent calculations along the Y dimension." +
        gpu_engaged_msg, tsm.getTestingStatus());
  check(synth_comz, RelationalOperator::EQUAL, synth_comz_ans, "Centers of mass computed using a "
        "MotionSweeper do not agree with basic, independent calculations along the Z dimension." +
        gpu_engaged_msg, tsm.getTestingStatus());
  check(synth_velx, RelationalOperator::EQUAL, synth_velx_ans, "Net velocities computed using a "
        "MotionSweeper do not agree with basic, independent calculations along the X dimension." +
        gpu_engaged_msg, tsm.getTestingStatus());
  check(synth_vely, RelationalOperator::EQUAL, synth_vely_ans, "Net velocities computed using a "
        "MotionSweeper do not agree with basic, independent calculations along the Y dimension." +
        gpu_engaged_msg, tsm.getTestingStatus());
  check(synth_velz, RelationalOperator::EQUAL, synth_velz_ans, "Net velocities computed using a "
        "MotionSweeper do not agree with basic, independent calculations along the Z dimension." +
        gpu_engaged_msg, tsm.getTestingStatus());

  // Move the systems to the origin and take away their net translational velocity
  removeCenterOfMassMotion(&mover_psw, mosw, gpu);
  std::vector<double> cenx(nsys, 0.0), ceny(nsys, 0.0), cenz(nsys, 0.0), zero_vec(nsys, 0.0);
  std::vector<double> netvx(nsys, 0.0), netvy(nsys, 0.0), netvz(nsys, 0.0);
  const TrajectoryKind apos = TrajectoryKind::POSITIONS;
  const TrajectoryKind velo = TrajectoryKind::VELOCITIES;
  for (int i = 0; i < nsys; i++) {
    const std::vector<double> centered_xyz = mover_ps.getInterlacedCoordinates(i, apos, insp_tier);
    const std::vector<double> velocity_xyz = mover_ps.getInterlacedCoordinates(i, velo, insp_tier);
    const size_t natom = host_mover_psw.atom_counts[i];
    const AtomGraph *i_ag = mover_ps.getSystemTopologyPointer(i);
    const std::vector<double> masses_i = i_ag->getAtomicMass<double>();
    for (size_t j = 0; j < natom; j++) {
      cenx[i] += centered_xyz[(3 * j)    ] * synth_masses[i][j];
      ceny[i] += centered_xyz[(3 * j) + 1] * synth_masses[i][j];
      cenz[i] += centered_xyz[(3 * j) + 2] * synth_masses[i][j];
      netvx[i] += velocity_xyz[(3 * j)    ] * synth_masses[i][j];
      netvy[i] += velocity_xyz[(3 * j) + 1] * synth_masses[i][j];
      netvz[i] += velocity_xyz[(3 * j) + 2] * synth_masses[i][j];
    }
    cenx[i] /= tmass[i];
    ceny[i] /= tmass[i];
    cenz[i] /= tmass[i];
    netvx[i] /= tmass[i];
    netvy[i] /= tmass[i];
    netvz[i] /= tmass[i];
  }
  check(cenx, RelationalOperator::EQUAL, zero_vec, "The MotionSweeper did not properly guide "
        "system centering along the X axis." + gpu_engaged_msg, tsm.getTestingStatus());
  check(ceny, RelationalOperator::EQUAL, zero_vec, "The MotionSweeper did not properly guide "
        "system centering along the Y axis." + gpu_engaged_msg, tsm.getTestingStatus());
  check(cenz, RelationalOperator::EQUAL, zero_vec, "The MotionSweeper did not properly guide "
        "system centering along the Z axis." + gpu_engaged_msg, tsm.getTestingStatus());
  check(netvx, RelationalOperator::EQUAL, zero_vec, "The MotionSweeper did not properly zero "
        "motion along the X axis." + gpu_engaged_msg, tsm.getTestingStatus());
  check(netvy, RelationalOperator::EQUAL, zero_vec, "The MotionSweeper did not properly zero "
        "motion along the Y axis." + gpu_engaged_msg, tsm.getTestingStatus());
  check(netvz, RelationalOperator::EQUAL, zero_vec, "The MotionSweeper did not properly zero "
        "motion along the Z axis." + gpu_engaged_msg, tsm.getTestingStatus());

  // Calculate the angular momentum and inertial moment
  accumulateAngularMomentum(&mosw, mover_auk, mover_psw, gpu);
  std::vector<double> angv_ans(3 * nsys, 0.0), angv_result(3 * nsys, 0.0);
  std::vector<double> itns_ans(9 * nsys, 0.0), itns_result(9 * nsys, 0.0);
  for (int i = 0; i < nsys; i++) {
    const std::vector<double> centered_xyz = mover_ps.getInterlacedCoordinates(i, apos, insp_tier);
    const std::vector<double> velocity_xyz = mover_ps.getInterlacedCoordinates(i, velo, insp_tier);
    const size_t natom = host_mover_psw.atom_counts[i];
    const AtomGraph *i_ag = mover_ps.getSystemTopologyPointer(i);
    const std::vector<double> masses_i = i_ag->getAtomicMass<double>();
    double r[3], v[3], rcv[3];
    double sum_ang[3] = { 0.0, 0.0, 0.0 };
    for (int j = 0; j < natom; j++) {
      r[0] = centered_xyz[(3 * j)    ];
      r[1] = centered_xyz[(3 * j) + 1];
      r[2] = centered_xyz[(3 * j) + 2];
      v[0] = velocity_xyz[(3 * j)    ];
      v[1] = velocity_xyz[(3 * j) + 1];
      v[2] = velocity_xyz[(3 * j) + 2];
      const double it_xx = r[0] * r[0] * masses_i[j];
      const double it_xy = r[0] * r[1] * masses_i[j];
      const double it_xz = r[0] * r[2] * masses_i[j];
      const double it_yy = r[1] * r[1] * masses_i[j];
      const double it_yz = r[1] * r[2] * masses_i[j];
      const double it_zz = r[2] * r[2] * masses_i[j];
      crossProduct(r, v, rcv);
      for (int k = 0; k < 3; k++) {
        sum_ang[k] += rcv[k] * masses_i[j];
      }
      itns_ans[(9 * i)    ] += it_yy + it_zz;
      itns_ans[(9 * i) + 1] -= it_xy;
      itns_ans[(9 * i) + 2] -= it_xz;
      itns_ans[(9 * i) + 3] -= it_xy;
      itns_ans[(9 * i) + 4] += it_xx + it_zz;
      itns_ans[(9 * i) + 5] -= it_yz;
      itns_ans[(9 * i) + 6] -= it_xz;
      itns_ans[(9 * i) + 7] -= it_yz;
      itns_ans[(9 * i) + 8] += it_xx + it_yy;
    }
    std::vector<double> itns_i(9), itns_inv(9);
    for (int j = 0; j < 9; j++) {
      itns_i[j] = itns_ans[(9 * i) + j];
    }
    invertSquareMatrix(itns_i, &itns_inv);
    for (int j = 0; j < 3; j++) {
      angv_ans[(3 * i) + j] = (itns_inv[j    ] * sum_ang[0]) + (itns_inv[j + 3] * sum_ang[1]) +
                              (itns_inv[j + 6] * sum_ang[2]);
    }
    const std::vector<double> inrt_i = mos.getInertialTensor(i, insp_tier);
    for (int j = 0; j < 9; j++) {
      itns_result[(9 * i) + j] = inrt_i[j];
    }
    const double3 angv_i = mos.getAngularVelocity(i, insp_tier);
    angv_result[(3 * i)    ] = angv_i.x;
    angv_result[(3 * i) + 1] = angv_i.y;
    angv_result[(3 * i) + 2] = angv_i.z;
  }
  check(itns_result, RelationalOperator::EQUAL, Approx(itns_ans).margin(1.7e-6), "Inertial "
        "tensors computed with a MotionSweeper do not match those computed using a basic "
        "approach." + gpu_engaged_msg, tsm.getTestingStatus());
  check(angv_result, RelationalOperator::EQUAL, angv_ans, "Net system angular velocities computed "
        "with a MotionSweeper do not match those computed using a basic approach." +
        gpu_engaged_msg, tsm.getTestingStatus());
  
  // Remove angular momentum, then check that the system is indeed no longer rotating and that it
  // still has no net translational velocity.
  removeAngularMomentum(&mover_psw, mosw, gpu);
  std::vector<double> residual_angv(3 * nsys, 0.0), residual_netv(3 * nsys, 0.0);
  for (int i = 0; i < nsys; i++) {
    const std::vector<double> centered_xyz = mover_ps.getInterlacedCoordinates(i, apos, insp_tier);
    const std::vector<double> velocity_xyz = mover_ps.getInterlacedCoordinates(i, velo, insp_tier);
    const size_t natom = host_mover_psw.atom_counts[i];
    const AtomGraph *i_ag = mover_ps.getSystemTopologyPointer(i);
    const std::vector<double> masses_i = i_ag->getAtomicMass<double>();
    double r[3], v[3], rcv[3];
    double sum_ang[3] = { 0.0, 0.0, 0.0 };
    for (int j = 0; j < natom; j++) {
      r[0] = centered_xyz[(3 * j)    ];
      r[1] = centered_xyz[(3 * j) + 1];
      r[2] = centered_xyz[(3 * j) + 2];
      v[0] = velocity_xyz[(3 * j)    ];
      v[1] = velocity_xyz[(3 * j) + 1];
      v[2] = velocity_xyz[(3 * j) + 2];
      crossProduct(r, v, rcv);
      for (int k = 0; k < 3; k++) {
        residual_netv[(3 * i) + k] += v[k] * masses_i[j];
        sum_ang[k] += rcv[k] * masses_i[j];
      }
    }
    for (int j = 0; j < 3; j++) {
      residual_netv[(3 * i) + j] /= tmass[i];
    }
    const std::vector<double> inrt_i = mos.getInertialTensor(i, insp_tier);
    std::vector<double> inv_inrt(9);
    invertSquareMatrix(inrt_i, &inv_inrt);
    for (int j = 0; j < 3; j++) {
      residual_angv[(3 * i) + j] = (inv_inrt[    j] * sum_ang[0]) +
                                   (inv_inrt[3 + j] * sum_ang[1]) +
                                   (inv_inrt[6 + j] * sum_ang[2]);
    }
  }
  check(residual_netv, RelationalOperator::EQUAL, std::vector<double>(3 * nsys, 0.0), "Removal of "
        "angular momentum leaves residual translational velocity." + gpu_engaged_msg,
        tsm.getTestingStatus());
  check(residual_angv, RelationalOperator::EQUAL, std::vector<double>(3 * nsys, 0.0), "Removal of "
        "angular momentum leaves residual angular velocity." + gpu_engaged_msg,
        tsm.getTestingStatus());
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Engage the testing environment
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
  StopWatch timer;

  // Engage the GPU
#ifdef STORMM_USE_HPC
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  const Hybrid<int> array_to_trigger_gpu_mapping(1);
#endif
  // Section 1: Test a harmonic spring with weights on either end.
  section("Weighted harmonic spring");

  // Section 2: Test momentum and net displacement removal
  section("Centering and momentum purge");

  // Section 3: Test the Thermostat class
  section("Thermostat unit tests");
  
  // Section 4: Test a molecule with all flexible bonds
  section("Flexible molecule");

  // Section 5: Test the Andersen thermostat
  section("Andersen thermostat");
  
  // Create the harmonic spring connecting two point masses, a 19.0 kcal/mol-A^2 spring with
  // equilibrium length connecting masses of 1.6 and 2.7 Daltons.  Simulate it for 1000 steps in
  // the isoenergetic ensemble.
  section(1);
  Dumbell d(1.6, 2.7, 10.0, 1.5);
  Xoshiro256ppGenerator xrs;
  std::vector<double3> init_r(2), init_v(2);
  for (int i = 0; i < 2; i++) {
    const double delta = static_cast<double>(i) * 2.0;
    init_r[i].x = xrs.gaussianRandomNumber() + delta;
    init_r[i].y = xrs.gaussianRandomNumber() + delta;
    init_r[i].z = xrs.gaussianRandomNumber() + delta;
    init_v[i].x = 0.05 * xrs.gaussianRandomNumber();
    init_v[i].y = 0.05 * xrs.gaussianRandomNumber();
    init_v[i].z = 0.05 * xrs.gaussianRandomNumber();
  }
  const double3 momentum = d.computeNetMomentum(init_v);
  const double tm = d.getFirstParticleMass() + d.getSecondParticleMass();
  const double3 subvel = { momentum.x / tm, momentum.y / tm, momentum.z / tm }; 
  for (int i = 0; i < 2; i++) {
    init_v[i].x -= subvel.x;
    init_v[i].y -= subvel.y;
    init_v[i].z -= subvel.z;
  }
  const int isonrg_sim_times = timer.addCategory("Iso-E dumbell simulations");
  timer.assignTime(0);
  double fs_two_var, fs_one_var, fs_hlf_var;
  const int nseg = 2000;
  const std::vector<double4> fs_two = isoenergetic(d, 2.0, (nseg + 1) * 1000, 1000, init_r, init_v,
                                                   &fs_two_var);
  const std::vector<double4> fs_one = isoenergetic(d, 1.0, (nseg + 1) * 2000, 2000, init_r, init_v,
                                                   &fs_one_var);
  const std::vector<double4> fs_hlf = isoenergetic(d, 0.5, (nseg + 1) * 4000, 4000, init_r, init_v,
                                                   &fs_hlf_var);
  timer.assignTime(isonrg_sim_times);
  const std::vector<double> fs_var = { fs_two_var, fs_one_var, fs_hlf_var };
  const std::vector<double> fs_var_ans = { 0.05585467, 0.01447308, 0.003689183 };
  check(fs_var, RelationalOperator::EQUAL, fs_var_ans, "The variance of total energies in "
        "simulations of a dumbell did not decrease as the square of the time step.  This could "
        "indicate that the integrator is not conserving energy.");

  // Create a handful of test systems
  const std::vector<std::string> mols = { "drug_example_iso", "bromobenzene_vs_iso", "hydrogenate",
                                          "drug_example_vs_iso", "med_1", "unsaturated_ring",
                                          "trpcage" };
  const char osc = osSeparator();
  const std::string base_top_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  TestSystemManager tsm(base_top_name, "top", mols, base_crd_name, "inpcrd", mols);

  // Create some simple systems and test momentum removal.
  section(2);
  runMomentumTests();  

  // Check momentum removal in the context of the synthesis.  Build up vectors of random
  // perturbations for each system, which can then be used to directly interrogate the recentering
  // and momentum removal.
  testMotionSweeper(tsm, &xrs);
#ifdef STORMM_USE_HPC
  // Make a new random number generator in order to repeat the tests on the GPU, so as to ensure
  // that subsequent tests using xrs will be consistent whether the GPU tests engage or not.
  Xoshiro256ppGenerator xrs_ii(45038194);
  testMotionSweeper(tsm, &xrs_ii, gpu);
#endif

  // Test the Thermostat class
  section(3);
  const AtomGraph ag_zero = tsm.exportAtomGraph(0);
  Thermostat tst_trial(ag_zero, ThermostatKind::LANGEVIN);
  CHECK_THROWS_SOFT(tst_trial.setLangevinCollisionFrequency(3.0), "A Langevin thermostat was set "
                    "with an exceedingly high collision frequency.", tsm.getTestingStatus());
  tst_trial.setLangevinCollisionFrequency(0.003);
  tst_trial.setTimeStep(4.0);
  tst_trial.setRandomCacheDepth(6);
  check(tst_trial.getTimeStep(), RelationalOperator::EQUAL, 4.0, "A Thermostat class object does "
        "not return the correct time step.", tsm.getTestingStatus());
  check(tst_trial.getRandomCacheDepth(), RelationalOperator::EQUAL, 6, "A Thermostat class object "
        "does not return the correct cache depth.", tsm.getTestingStatus());
  check(tst_trial.getAtomCount(), RelationalOperator::EQUAL, ag_zero.getAtomCount(),
        "A Thermostat class object does not record the correct number of atoms.",
        tsm.getTestingStatus());
  tst_trial.setCacheConfiguration(PrecisionModel::DOUBLE);
  check(tst_trial.getCacheConfiguration() == PrecisionModel::DOUBLE, "The random number cache of "
        "a thermostat was not altered as requested.", tsm.getTestingStatus());
  tst_trial.initializeRandomStates(561249238, 25, HybridTargetLevel::HOST);
  PhaseSpace ps_zero = tsm.exportPhaseSpace(0);
  DynamicsControls tmp_dyncon;
  tmp_dyncon.setGeometricConstraints(ApplyConstraints::NO);
  velocityKickStart(&ps_zero, ag_zero, &tst_trial, tmp_dyncon, EnforceExactTemperature::NO);
  ScoreCard sc_zero(1, 1, 32);
  evalKineticEnergy(ps_zero, CoordinateCycle::BLACK, &sc_zero, ag_zero);
  computeTemperature(&sc_zero, ag_zero, false);
  check(sc_zero.reportInstantaneousStates(StateVariable::TEMPERATURE_ALL, 0),
        RelationalOperator::EQUAL, Approx(302.75093816).margin(1.0e-5), "The temperature of a "
        "small system kick-started with an Andersen thermostat, no geometry constraints, does "
        "not meet expectations.");

  // Create a few molecules in isolated boundary conditions and test their dynamics using the
  // velocity-Verlet integrator.
  section(4);
  std::vector<double> temperatures;
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    PhaseSpace ps = tsm.exportPhaseSpace(i);
    AtomGraph ag = tsm.exportAtomGraph(i);
    double nrg_cons_tol, eqi_part_tol;
    if (mols[i] == "hydrogenate") {
      nrg_cons_tol = 1.6e-2;
      eqi_part_tol = 3.0e-2;
    }
    else if (mols[i] == "trpcage") {
      nrg_cons_tol = 1.8e-1;
      eqi_part_tol = 1.4e-1;
    }
    else {
      eqi_part_tol = 3.0e-2;
      nrg_cons_tol = 4.0e-3;
    }
    const StaticExclusionMask se(ag);
    if (mols[i] == "hydrogenate" || mols[i] == "unsaturated_ring") {
      PhaseSpaceWriter psw = ps.data();
      addRandomNoise(&xrs, psw.xcrd, psw.ycrd, psw.zcrd, psw.natom, 0.005);
    }
    timer.assignTime(0);
    const int tmcat = timer.addCategory("MD on " + mols[i]);
    const double ta = mdTrial(&ps, ag, 8000, 0.25, 0.0, RattleMethod::SEQUENTIAL,
                              tsm.getTestingStatus(), 2000, eqi_part_tol, nrg_cons_tol);
    const double tb = mdTrial(&ps, ag, 8000, 0.50, 0.0, RattleMethod::SEQUENTIAL,
                              tsm.getTestingStatus(), 200, eqi_part_tol, nrg_cons_tol);
    const double tc = mdTrial(&ps, ag, 8000, 1.00, 0.0, RattleMethod::SEQUENTIAL,
                              tsm.getTestingStatus(), 2000, eqi_part_tol, nrg_cons_tol);
    const double td = mdTrial(&ps, ag, 8000, 2.00, 1.0e-7, RattleMethod::SEQUENTIAL,
                              tsm.getTestingStatus(), 2000, eqi_part_tol, nrg_cons_tol);
    const double te = mdTrial(&ps, ag, 8000, 2.00, 1.0e-7, RattleMethod::CENTER_SUM,
                              tsm.getTestingStatus(), 2000, eqi_part_tol, nrg_cons_tol);
    temperatures.push_back(ta);
    temperatures.push_back(tb);
    temperatures.push_back(tc);
    temperatures.push_back(td);
    temperatures.push_back(te);
    timer.assignTime(tmcat);
  }

  // Check the temperatures against established mean values
  const std::string ref_file = oe.getStormmSourcePath() + osc + "test" + osc +
                               "MolecularMechanics" + osc + "mdyn_output.m";
  const TestPriority do_snps = (getDrivePathType(ref_file) == DrivePathType::FILE) ?
                               TestPriority::CRITICAL : TestPriority::ABORT;
  snapshot(ref_file, polyNumericVector(temperatures), "vac_temp", 4.0, "Mean temperatures from "
           "short MD simulations of very small systems (conducted in vacuum, with "
           "double-precision arithmetic throughout) do not meet expectations.", oe.takeSnapshot(),
           1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::OVERWRITE, do_snps);
  std::vector<double> ref_t;
  switch (do_snps) {
  case TestPriority::CRITICAL:
    {
      const TextFile ref_data(ref_file);
      ref_t = doubleFromPolyNumeric(readSnapshot(ref_data, "vac_temp"));
    }
    break;
  case TestPriority::NON_CRITICAL:
  case TestPriority::ABORT:
    break;
  }
  check(pearson(ref_t, temperatures), RelationalOperator::EQUAL, Approx(1.0).margin(2.0e-2),
        "Mean temperatures obtained from a series of short MD simulations do not correlate with "
        "the expected values.  Deviations from platform to platform are expected, but the mean "
        "temperatures, reflecting the overall amounts of energy injected into each system by "
        "consistent initial conditions, should align.", do_snps);

  // Test the particle temperature initialization with the Andersen themostat.  Include restraints
  // and implicit solvent models.
  MinimizeControls mincon;
  mincon.setTotalCycles(250);
  mincon.setClashDampingCycles(50);
  mincon.setDiagnosticPrintFrequency(25);
  NeckGeneralizedBornTable ngb_tab;
  DynamicsControls dyncon;
  dyncon.setGeometricConstraints();
  dyncon.setCpuRattleMethod("center_sum");
  dyncon.setRattleTolerance(1.0e-7);
  dyncon.setDiagnosticPrintFrequency(20);
  dyncon.setStepCount(5440);
  const int nsys = tsm.getSystemCount();
  std::vector<double> bond_emin_delta(nsys), angl_emin_delta(nsys), dihe_emin_delta(nsys);
  std::vector<double> impr_emin_delta(nsys), elec_emin_delta(nsys), vdwa_emin_delta(nsys);
  std::vector<double> qq14_emin_delta(nsys), lj14_emin_delta(nsys), rstr_emin_delta(nsys);
  std::vector<double> gbrn_emin_delta(nsys), totl_emin_delta(nsys);
  std::vector<double> bond_md_ave(nsys), angl_md_ave(nsys), dihe_md_ave(nsys), impr_md_ave(nsys);
  std::vector<double> elec_md_ave(nsys), vdwa_md_ave(nsys), qq14_md_ave(nsys), lj14_md_ave(nsys);
  std::vector<double> rstr_md_ave(nsys), gbrn_md_ave(nsys), kine_md_ave(nsys), totl_md_ave(nsys);
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    PhaseSpace ps = tsm.exportPhaseSpace(i);
    AtomGraph ag = tsm.exportAtomGraph(i);
    const StaticExclusionMask se(ag);
      
    // If there are no virtual sites, add an implicit solvent model (GB and virtual sites are
    // seldom tested together in computational chemistry).
    if (ag.getVirtualSiteCount() == 0) {
      const ImplicitSolventModel ismod = static_cast<ImplicitSolventModel>(i % 6);
      AtomicRadiusSet rads;
      switch (ismod) {
      case ImplicitSolventModel::NONE:
        rads = AtomicRadiusSet::NONE;
        break;
      case ImplicitSolventModel::HCT_GB:
      case ImplicitSolventModel::OBC_GB:
      case ImplicitSolventModel::OBC_GB_II:
        rads = AtomicRadiusSet::MBONDI2;
        break;
      case ImplicitSolventModel::NECK_GB:
      case ImplicitSolventModel::NECK_GB_II:
        rads = AtomicRadiusSet::MBONDI3;
        break;
      }
      ag.setImplicitSolventModel(ismod, 80.0, 0.0, rads, ExceptionResponse::WARN);
    }

    // Include restraints of different orders between random atoms in the structure.
    RestraintApparatus ra(&ag);
    PhaseSpaceWriter psw = ps.data();
    const CoordinateFrameReader cfr(ps);
    std::vector<BoundedRestraint> br_adds;
    for (int j = 0; j < 4; j++) {

      // Add a position restraint
      const int atom_a = static_cast<double>(psw.natom) * xrs.uniformRandomNumber();
      br_adds.emplace_back(atom_a, &ag, cfr, 1.0, 2.0, 0.0, 0.5, 2.5, 8.5, atom_a);

      // Add an angle restraint
      const int atom_b = static_cast<double>(psw.natom / 4) * xrs.uniformRandomNumber();
      const int atom_c = atom_b + 3;
      const int atom_d = atom_c + 4;
      if (atom_d < psw.natom) {
        const double current_angl = angle(atom_b, atom_c, atom_d, ps);
        br_adds.emplace_back(atom_b, atom_c, atom_d, &ag, 4.0, 4.0, current_angl - 0.6,
                             current_angl - 0.1, current_angl - 0.1, current_angl + 0.4);
      }
    }
    ra.addRestraints(br_adds);
    ScoreCard sc = minimize(&ps, ag, ngb_tab, ra, se, mincon);

    // Check the energy minimization (for additional coverage of minimization protocols)
    sc.computeTotalEnergy();
    const std::vector<double> totl_m = sc.reportHistory(StateVariable::TOTAL_ENERGY, 0);
    const std::vector<double> bond_m = sc.reportHistory(StateVariable::BOND, 0);
    const std::vector<double> angl_m = sc.reportHistory(StateVariable::ANGLE, 0);
    const std::vector<double> dihe_m = sc.reportHistory(StateVariable::PROPER_DIHEDRAL, 0);
    const std::vector<double> impr_m = sc.reportHistory(StateVariable::IMPROPER_DIHEDRAL, 0);
    const std::vector<double> elec_m = sc.reportHistory(StateVariable::ELECTROSTATIC, 0);
    const std::vector<double> vdwa_m = sc.reportHistory(StateVariable::VDW, 0);
    const std::vector<double> qq14_m = sc.reportHistory(StateVariable::ELEC_ONE_FOUR, 0);
    const std::vector<double> lj14_m = sc.reportHistory(StateVariable::VDW_ONE_FOUR, 0);
    const std::vector<double> rstr_m = sc.reportHistory(StateVariable::RESTRAINT, 0);
    const std::vector<double> gbrn_m = sc.reportHistory(StateVariable::GENERALIZED_BORN, 0);
    const std::vector<double> kine_m = sc.reportHistory(StateVariable::KINETIC, 0);
    bond_emin_delta[i] = bond_m.back() - bond_m.front();
    angl_emin_delta[i] = angl_m.back() - angl_m.front();
    dihe_emin_delta[i] = dihe_m.back() - dihe_m.front();
    impr_emin_delta[i] = impr_m.back() - impr_m.front();
    elec_emin_delta[i] = elec_m.back() - elec_m.front();
    vdwa_emin_delta[i] = vdwa_m.back() - vdwa_m.front();
    qq14_emin_delta[i] = qq14_m.back() - qq14_m.front();
    lj14_emin_delta[i] = lj14_m.back() - lj14_m.front();
    rstr_emin_delta[i] = rstr_m.back() - rstr_m.front();
    gbrn_emin_delta[i] = gbrn_m.back() - gbrn_m.front();
    totl_emin_delta[i] = totl_m.back() - totl_m.front();
    
    // Copy the minimized coordinates into the alternate stage of the coordinate cycle.
    for (int i = 0; i < psw.natom; i++) {
      psw.xalt[i] = psw.xcrd[i];
      psw.yalt[i] = psw.ycrd[i];
      psw.zalt[i] = psw.zcrd[i];
    }
    
    // Submit the minimized structure to be kick-started into dynamics with initialization by an
    // Andersen thermostat cycle.
    Thermostat tst(ag, ThermostatKind::LANGEVIN, 300.0);
    tst.setLangevinCollisionFrequency(0.003);
    tst.setTimeStep(1.0);
    tst.setGeometryConstraints(dyncon.constrainGeometry());
    velocityKickStart(&ps, ag, &tst, dyncon);
    
    // Check that the temperature is now set to exactly 300 Kelvin
    double init_temp;
    std::string has_cnst;
    switch (dyncon.constrainGeometry()) {
    case ApplyConstraints::YES:
      init_temp = computeTemperature<double>(ps, ag, ApplyConstraints::YES);
      has_cnst = "with";
      break;
    case ApplyConstraints::NO:
      init_temp = computeTemperature<double>(ps, ag, ApplyConstraints::NO);
      has_cnst = "without";
      break;
    }
    check(init_temp, RelationalOperator::EQUAL, Approx(300.0).margin(1.0e-4), "The temperature of "
          "a system kick-started " + has_cnst + " constraints does not meet expectations.",
          tsm.getTestingStatus());

    // Run dynamics at constant energy
    ScoreCard scmd(1, 101, 32);
    dynamics(&ps, &tst, &scmd, ag, ngb_tab, se, ra, dyncon, 0);
    scmd.computeTotalEnergy();
    const std::vector<double> totl_e = scmd.reportHistory(StateVariable::TOTAL_ENERGY, 0);
    const std::vector<double> bond_e = scmd.reportHistory(StateVariable::BOND, 0);
    const std::vector<double> angl_e = scmd.reportHistory(StateVariable::ANGLE, 0);
    const std::vector<double> dihe_e = scmd.reportHistory(StateVariable::PROPER_DIHEDRAL, 0);
    const std::vector<double> impr_e = scmd.reportHistory(StateVariable::IMPROPER_DIHEDRAL, 0);
    const std::vector<double> elec_e = scmd.reportHistory(StateVariable::ELECTROSTATIC, 0);
    const std::vector<double> vdwa_e = scmd.reportHistory(StateVariable::VDW, 0);
    const std::vector<double> qq14_e = scmd.reportHistory(StateVariable::ELEC_ONE_FOUR, 0);
    const std::vector<double> lj14_e = scmd.reportHistory(StateVariable::VDW_ONE_FOUR, 0);
    const std::vector<double> rstr_e = scmd.reportHistory(StateVariable::RESTRAINT, 0);
    const std::vector<double> gbrn_e = scmd.reportHistory(StateVariable::GENERALIZED_BORN, 0);
    const std::vector<double> kine_e = scmd.reportHistory(StateVariable::KINETIC, 0);
    bond_md_ave[i] = mean(&bond_e[16], 256);
    angl_md_ave[i] = mean(&angl_e[16], 256);
    dihe_md_ave[i] = mean(&dihe_e[16], 256);
    impr_md_ave[i] = mean(&impr_e[16], 256);
    elec_md_ave[i] = mean(&elec_e[16], 256);
    vdwa_md_ave[i] = mean(&vdwa_e[16], 256);
    qq14_md_ave[i] = mean(&qq14_e[16], 256);
    lj14_md_ave[i] = mean(&lj14_e[16], 256);
    rstr_md_ave[i] = mean(&rstr_e[16], 256);
    gbrn_md_ave[i] = mean(&gbrn_e[16], 256);
    kine_md_ave[i] = mean(&kine_e[16], 256);
    totl_md_ave[i] = mean(&totl_e[16], 256);
  }
  const TestPriority mean_nrg_comparisons = (do_snps == TestPriority::CRITICAL) ?
                                            TestPriority::NON_CRITICAL : TestPriority::ABORT;
  snapshot(ref_file, polyNumericVector(totl_emin_delta), "totl_min", 4.0, "Mean energy "
           "differences energy minimizations of very small systems (conducted in various implicit "
           "solvent models, with double-precision arithmetic throughout) do not meet "
           "expectations.", oe.takeSnapshot(), 1.0e-4, NumberFormat::STANDARD_REAL,
           PrintSituation::APPEND, mean_nrg_comparisons);
  snapshot(ref_file, polyNumericVector(totl_md_ave), "totl_ave", 4.0, "Mean energies obtained "
           "from short simulations of very small systems (conducted in various implicit "
           "solvent models, with double-precision arithmetic throughout) do not meet "
           "expectations.", oe.takeSnapshot(), 2.0e-3, NumberFormat::STANDARD_REAL,
           PrintSituation::APPEND, mean_nrg_comparisons);
  std::vector<double> ref_totl_emin_delta, ref_totl_md_ave;
  switch (do_snps) {
  case TestPriority::CRITICAL:
    {
      const TextFile ref_data(ref_file);
      ref_totl_emin_delta = doubleFromPolyNumeric(readSnapshot(ref_data, "totl_min"));
      ref_totl_md_ave = doubleFromPolyNumeric(readSnapshot(ref_data, "totl_ave"));
    }
    break;
  case TestPriority::NON_CRITICAL:
  case TestPriority::ABORT:
    break;
  }

  // Check the energy minimization deltas, and the MD energy averages
  const TestPriority emin_comparisons = (do_snps == TestPriority::CRITICAL) ?
                                        TestPriority::NON_CRITICAL : TestPriority::ABORT;
  snapshot(ref_file, polyNumericVector(bond_emin_delta), "bond_min", 1.0e-1, "Bond energy deltas "
           "obtained from minimizations do not meet expectations.", oe.takeSnapshot(), 1.0e-4,
           NumberFormat::STANDARD_REAL, PrintSituation::APPEND, emin_comparisons);
  snapshot(ref_file, polyNumericVector(angl_emin_delta), "angl_min", 1.0e-1, "Angle energy deltas "
           "obtained from minimizations do not meet expectations.", oe.takeSnapshot(), 1.0e-4,
           NumberFormat::STANDARD_REAL, PrintSituation::APPEND, emin_comparisons);
  snapshot(ref_file, polyNumericVector(dihe_emin_delta), "dihe_min", 1.0e-1, "Dihedral energy "
           "deltas obtained from minimizations do not meet expectations.", oe.takeSnapshot(),
           1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, emin_comparisons);
  snapshot(ref_file, polyNumericVector(impr_emin_delta), "impr_min", 1.0e-1, "Improper energy "
           "deltas obtained from minimizations do not meet expectations.", oe.takeSnapshot(),
           1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, emin_comparisons);
  snapshot(ref_file, polyNumericVector(elec_emin_delta), "elec_min", 1.0e-1, "Electrostatic "
           "energy deltas obtained from minimizations do not meet expectations.",
           oe.takeSnapshot(), 1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND,
           emin_comparisons);
  snapshot(ref_file, polyNumericVector(vdwa_emin_delta), "vdwa_min", 1.0e-1, "Van-der Waals "
           "energy deltas obtained from minimizations do not meet expectations.",
           oe.takeSnapshot(), 1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND,
           emin_comparisons);
  snapshot(ref_file, polyNumericVector(qq14_emin_delta), "qq14_min", 1.0e-1, "Electrostatic 1-4 "
           "energy deltas obtained from minimizations do not meet expectations.",
           oe.takeSnapshot(), 1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND,
           emin_comparisons);
  snapshot(ref_file, polyNumericVector(lj14_emin_delta), "lj14_min", 1.0e-1, "Lennard-Jones 1-4 "
           "energy deltas obtained from minimizations do not meet expectations.",
           oe.takeSnapshot(), 1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND,
           emin_comparisons);
  snapshot(ref_file, polyNumericVector(rstr_emin_delta), "rstr_min", 1.0e-1, "Restraint energy "
           "deltas obtained from minimizations do not meet expectations.", oe.takeSnapshot(),
           1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, emin_comparisons);
  snapshot(ref_file, polyNumericVector(gbrn_emin_delta), "gbrn_min", 1.0e-1, "Generalized Born "
           "energy deltas obtained from minimizations do not meet expectations.",
           oe.takeSnapshot(), 1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND,
           emin_comparisons);
  snapshot(ref_file, polyNumericVector(totl_emin_delta), "totl_min", 1.0e-1, "Total energy "
           "deltas obtained from minimizations do not meet expectations.", oe.takeSnapshot(),
           1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, emin_comparisons);
  check(pearson(totl_emin_delta, ref_totl_emin_delta), RelationalOperator::EQUAL,
        Approx(1.0).margin(1.0e-2), "The total energies released by relaxation procedures applied "
        "to a variety of small systems should align with reference values, even if the exact "
        "pathways taken by the energy minimization and the local minima they fall into contain "
        "slight differences depending on the platform.", do_snps);
  
  // Check the energy minimization deltas, and the MD energy averages
  const TestPriority mdrun_comparisons = (do_snps == TestPriority::CRITICAL) ?
                                         TestPriority::NON_CRITICAL : TestPriority::ABORT;
  snapshot(ref_file, polyNumericVector(bond_md_ave), "bond_ave", 1.0e-1, "Bond energy averages "
           "obtained from MD simulations do not meet expectations.", oe.takeSnapshot(), 1.0e-4,
           NumberFormat::STANDARD_REAL, PrintSituation::APPEND, mdrun_comparisons);
  snapshot(ref_file, polyNumericVector(angl_md_ave), "angl_ave", 1.0e-1, "Angle energy averages "
           "obtained from MD simulations do not meet expectations.", oe.takeSnapshot(), 1.0e-4,
           NumberFormat::STANDARD_REAL, PrintSituation::APPEND, mdrun_comparisons);
  snapshot(ref_file, polyNumericVector(dihe_md_ave), "dihe_ave", 1.0e-1, "Dihedral energy "
           "averages obtained from MD simulations do not meet expectations.", oe.takeSnapshot(),
           1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, mdrun_comparisons);
  snapshot(ref_file, polyNumericVector(impr_md_ave), "impr_ave", 1.0e-1, "Improper energy "
           "averages obtained from MD simulations do not meet expectations.", oe.takeSnapshot(),
           1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, mdrun_comparisons);
  snapshot(ref_file, polyNumericVector(elec_md_ave), "elec_ave", 1.0e-1, "Electrostatic "
           "energy averages obtained from MD simulations do not meet expectations.",
           oe.takeSnapshot(), 1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND,
           mdrun_comparisons);
  snapshot(ref_file, polyNumericVector(vdwa_md_ave), "vdwa_ave", 1.0e-1, "Van-der Waals "
           "energy averages obtained from MD simulations do not meet expectations.",
           oe.takeSnapshot(), 1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND,
           mdrun_comparisons);
  snapshot(ref_file, polyNumericVector(qq14_md_ave), "qq14_ave", 1.0e-1, "Electrostatic 1-4 "
           "energy averages obtained from MD simulations do not meet expectations.",
           oe.takeSnapshot(), 1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND,
           mdrun_comparisons);
  snapshot(ref_file, polyNumericVector(lj14_md_ave), "lj14_ave", 1.0e-1, "Lennard-Jones 1-4 "
           "energy averages obtained from MD simulations do not meet expectations.",
           oe.takeSnapshot(), 1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND,
           mdrun_comparisons);
  snapshot(ref_file, polyNumericVector(rstr_md_ave), "rstr_ave", 1.0e-1, "Restraint energy "
           "averages obtained from MD simulations do not meet expectations.", oe.takeSnapshot(),
           1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, mdrun_comparisons);
  snapshot(ref_file, polyNumericVector(gbrn_md_ave), "gbrn_ave", 1.0e-1, "Generalized Born "
           "energy averages obtained from MD simulations do not meet expectations.",
           oe.takeSnapshot(), 1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND,
           mdrun_comparisons);
  snapshot(ref_file, polyNumericVector(gbrn_md_ave), "gbrn_ave", 1.0e-1, "Kinetic energy "
           "averages obtained from MD simulations do not meet expectations.", oe.takeSnapshot(),
           1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, mdrun_comparisons);
  snapshot(ref_file, polyNumericVector(totl_md_ave), "totl_ave", 1.5e-1, "Total energy "
           "averages obtained from MD simulations do not meet expectations.", oe.takeSnapshot(),
           1.0e-4, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, mdrun_comparisons);
  check(pearson(totl_md_ave, ref_totl_md_ave), RelationalOperator::EQUAL,
        Approx(1.0).margin(1.6e-2), "The total energies maintained by conservative MD protocols "
        "applied to a variety of small systems should align with reference values.", do_snps);

  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
