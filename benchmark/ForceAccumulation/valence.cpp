#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/fixed_precision.h"
#include "../../src/Constants/hpc_bounds.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/Trajectory/trajectory_enumerators.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::constants::ExceptionResponse;
using stormm::numerics::default_globalpos_scale_lf;
using stormm::numerics::default_globalpos_scale_f;
using stormm::numerics::default_inverse_globalpos_scale_lf;
using stormm::numerics::default_inverse_globalpos_scale_f;
using stormm::data_types::double_type_index;
#ifndef STORMM_USE_HPC
using stormm::data_types::double3;
using stormm::data_types::float2;
using stormm::data_types::float3;
#endif
using stormm::data_types::llint;
using stormm::data_types::llint3;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::energy::EvaluateForce;
using stormm::energy::ScoreCard;
using stormm::errors::terminalFormat;
using stormm::parse::PolyNumeric;
using stormm::random::Xoshiro256ppGenerator;
using stormm::topology::AtomGraph;
using stormm::topology::ValenceKit;
using stormm::trajectory::CoordinateFileKind;
using stormm::trajectory::PhaseSpace;
using stormm::trajectory::TrajectoryKind;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Compute the force between two particles given a set of harmonic bond parameters.
//
// Arguments:
//   crd1:   Coordinates of the first particle
//   crd2:   Coordinates of the second particle
//   equil:  Bond equilibrium length 
//   stiff:  Bond stiffness constant
//-------------------------------------------------------------------------------------------------
template <typename T, typename T3>
T3 bond_force(T3 crd1, T3 crd2, T equil, T stiff) {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const T dx = crd2.x - crd1.x;
  const T dy = crd2.y - crd1.y;
  const T dz = crd2.z - crd1.z;
  const T r = (ct == double_type_index) ? sqrt((dx * dx) + (dy * dy) + (dz * dz)) :
                                          sqrtf((dx * dx) + (dy * dy) + (dz * dz));

  // Select a random equilibrium length, based on the target length, and a stiffness constant
  const T dl = r - equil;
  const T fmag = 2.0 * stiff * dl / r;
  T3 result;
  result.x = fmag * dx;
  result.y = fmag * dy;
  result.z = fmag * dz;
  return result;
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Baseline variables
  TestEnvironment oe(argc, argv);

  // Perform arithmetic operations similar to a single bond force computation for a range of
  // distances and stiffness constants that emulates what would be encountered in a typical
  // simulation.  Compute the quantities with coordinates and parameter constants in double
  // precision, then in single precision, then in split single precision with fixed-precision
  // intermediaries and accumulation.
  const int nsample = 256;
  const int ndistance = 1200;
  const int npts = nsample * ndistance;
  const double distance_discretization = 0.001;
  Xoshiro256ppGenerator xsr_rng(90384011);
  std::vector<double3> crd1(npts), crd2(npts), frc(npts), frc_with_fequil(npts);
  std::vector<float3> f_crd1(npts), f_crd2(npts), f_frc(npts), flli_frc(npts), f2lli_frc(npts);
  std::vector<llint3> lli_crd1(npts), lli_crd2(npts);
  std::vector<double> equil(npts), stiff(npts);
  std::vector<float> f_equil(npts), f_stiff(npts);
  std::vector<float2> f2_equil(npts), f2_stiff(npts);
  for (int i = 0; i < npts; i++) {

    // Select random coordinates
    crd1[i].x = xsr_rng.gaussianRandomNumber();
    crd1[i].y = xsr_rng.gaussianRandomNumber();
    crd1[i].z = xsr_rng.gaussianRandomNumber();
    crd2[i].x = xsr_rng.gaussianRandomNumber();
    crd2[i].y = xsr_rng.gaussianRandomNumber();
    crd2[i].z = xsr_rng.gaussianRandomNumber();

    // Extend or contract the distance along the same direction to meet a particular value
    const double target_distance = 0.8 + ((static_cast<double>(i / nsample) +
                                           xsr_rng.uniformRandomNumber()) *
                                          distance_discretization);
    const double pdx = crd2[i].x - crd1[i].x;
    const double pdy = crd2[i].y - crd1[i].y;
    const double pdz = crd2[i].z - crd1[i].z;
    const double factor = target_distance / sqrt((pdx * pdx) + (pdy * pdy) + (pdz * pdz));
    crd1[i].x += (xsr_rng.uniformRandomNumber() - 0.5) * 80.0;
    crd1[i].y += (xsr_rng.uniformRandomNumber() - 0.5) * 80.0;
    crd1[i].z += (xsr_rng.uniformRandomNumber() - 0.5) * 80.0;
    crd2[i].x = crd1[i].x + (factor * pdx);
    crd2[i].y = crd1[i].y + (factor * pdy);
    crd2[i].z = crd1[i].z + (factor * pdz);
    stiff[i] = 350.0 + (250.0 * xsr_rng.uniformRandomNumber());
    equil[i] = target_distance + ((50.0 / stiff[i]) * xsr_rng.gaussianRandomNumber());
    
    // Compute the force in double precision
    frc[i] = bond_force(crd1[i], crd2[i], equil[i], stiff[i]);
    
    // Perform the computation with all coordinates reduced to floating-point numbers
    f_crd1[i].x = crd1[i].x;
    f_crd1[i].y = crd1[i].y;
    f_crd1[i].z = crd1[i].z;
    f_crd2[i].x = crd2[i].x;
    f_crd2[i].y = crd2[i].y;
    f_crd2[i].z = crd2[i].z;
    f_equil[i] = equil[i];
    f_stiff[i] = stiff[i];    
    f_frc[i] = bond_force(f_crd1[i], f_crd2[i], f_equil[i], f_stiff[i]);
    frc_with_fequil[i] = bond_force(crd1[i], crd2[i], static_cast<double>(f_equil[i]), stiff[i]);

    // Perform the computation with all coordinates reduced to long long integers
    lli_crd1[i].x = crd1[i].x * default_globalpos_scale_lf;
    lli_crd1[i].y = crd1[i].y * default_globalpos_scale_lf;
    lli_crd1[i].z = crd1[i].z * default_globalpos_scale_lf;
    lli_crd2[i].x = crd2[i].x * default_globalpos_scale_lf;
    lli_crd2[i].y = crd2[i].y * default_globalpos_scale_lf;
    lli_crd2[i].z = crd2[i].z * default_globalpos_scale_lf;
    {
      const float dx = static_cast<float>(lli_crd2[i].x - lli_crd1[i].x) *
                       default_inverse_globalpos_scale_f;
      const float dy = static_cast<float>(lli_crd2[i].y - lli_crd1[i].y) *
                       default_inverse_globalpos_scale_f;
      const float dz = static_cast<float>(lli_crd2[i].z - lli_crd1[i].z) *
                       default_inverse_globalpos_scale_f;
      const float r = sqrtf((dx * dx) + (dy * dy) + (dz * dz));
      const float dl = r - f_equil[i];
      const float fmag = 2.0 * f_stiff[i] * dl / r;
      flli_frc[i].x = fmag * dx;
      flli_frc[i].y = fmag * dy;
      flli_frc[i].z = fmag * dz;
    }

    // Perform the calculation with long long integer arithmetic up to the sqrt
    {
      const llint idx = lli_crd2[i].x - lli_crd1[i].x;
      const llint idy = lli_crd2[i].y - lli_crd1[i].y;
      const llint idz = lli_crd2[i].z - lli_crd1[i].z;
      if (llabs(idx) < 1600000000LL && llabs(idy) < 1600000000LL && llabs(idz) < 1600000000LL) {
        const llint ir2 = (idx * idx) + (idy * idy) + (idz * idz);
        const double r = sqrt(static_cast<double>(ir2)) * default_inverse_globalpos_scale_lf;
        const float dl = static_cast<float>(r - equil[i]);
        const float fmag = 2.0 * f_stiff[i] * dl / r;
        const float dx = static_cast<float>(idx) * default_inverse_globalpos_scale_f;
        const float dy = static_cast<float>(idy) * default_inverse_globalpos_scale_f;
        const float dz = static_cast<float>(idz) * default_inverse_globalpos_scale_f;
        f2lli_frc[i].x = fmag * dx;
        f2lli_frc[i].y = fmag * dy;
        f2lli_frc[i].z = fmag * dz;
      }
      else {
        const float dx = static_cast<float>(lli_crd2[i].x - lli_crd1[i].x) *
                         default_inverse_globalpos_scale_f;
        const float dy = static_cast<float>(lli_crd2[i].y - lli_crd1[i].y) *
                         default_inverse_globalpos_scale_f;
        const float dz = static_cast<float>(lli_crd2[i].z - lli_crd1[i].z) *
                         default_inverse_globalpos_scale_f;
        const float r = sqrtf((dx * dx) + (dy * dy) + (dz * dz));
        const float dl = r - f_equil[i];
        const float fmag = 2.0 * f_stiff[i] * dl / r;
        f2lli_frc[i].x = fmag * dx;
        f2lli_frc[i].y = fmag * dy;
        f2lli_frc[i].z = fmag * dz;
      }
    }
  }
  
  // Compute and report statistics from the initial experiment
  double mue_f            = 0.0;
  double mue_f_mag        = 0.0;
  double mue_flli         = 0.0;
  double mue_flli_mag     = 0.0;
  double mue_f2lli        = 0.0;
  double mue_f2lli_mag    = 0.0;
  for (int i = 0; i < npts; i++) {
    const double dx = static_cast<double>(f_frc[i].x) - frc[i].x;
    const double dy = static_cast<double>(f_frc[i].y) - frc[i].y;
    const double dz = static_cast<double>(f_frc[i].z) - frc[i].z;
    mue_f_mag += sqrt((dx * dx) + (dy * dy) + (dz * dz));
    mue_f += fabs(static_cast<double>(f_frc[i].x) - frc[i].x);
    mue_f += fabs(static_cast<double>(f_frc[i].y) - frc[i].y);
    mue_f += fabs(static_cast<double>(f_frc[i].z) - frc[i].z);
  }
  for (int i = 0; i < npts; i++) {
    const double dx = static_cast<double>(flli_frc[i].x) - frc[i].x;
    const double dy = static_cast<double>(flli_frc[i].y) - frc[i].y;
    const double dz = static_cast<double>(flli_frc[i].z) - frc[i].z;
    mue_flli_mag += sqrt((dx * dx) + (dy * dy) + (dz * dz));
    mue_flli += fabs(static_cast<double>(flli_frc[i].x) - frc[i].x);
    mue_flli += fabs(static_cast<double>(flli_frc[i].y) - frc[i].y);
    mue_flli += fabs(static_cast<double>(flli_frc[i].z) - frc[i].z);
  }
  for (int i = 0; i < npts; i++) {
    const double dx = static_cast<double>(f2lli_frc[i].x) - frc[i].x;
    const double dy = static_cast<double>(f2lli_frc[i].y) - frc[i].y;
    const double dz = static_cast<double>(f2lli_frc[i].z) - frc[i].z;
    mue_f2lli_mag += sqrt((dx * dx) + (dy * dy) + (dz * dz));
    mue_f2lli += fabs(static_cast<double>(f2lli_frc[i].x) - frc[i].x);
    mue_f2lli += fabs(static_cast<double>(f2lli_frc[i].y) - frc[i].y);
    mue_f2lli += fabs(static_cast<double>(f2lli_frc[i].z) - frc[i].z);
  }
  const double dnpts = static_cast<double>(npts);
  mue_f            /= 3.0 * dnpts;
  mue_flli         /= 3.0 * dnpts;
  mue_f2lli        /= 3.0 * dnpts;
  mue_f_mag        /= dnpts;
  mue_flli_mag     /= dnpts;
  mue_f2lli_mag    /= dnpts;
  std::string output;
  output = terminalFormat("Bond parameters and coordinates represented in pure 32-bit floating "
                          "point numbers, all calculations done in 32-bit floating point "
                          "arithmetic:", nullptr, nullptr, 0, 1, 1);
  printf("%s\n", output.c_str());
  printf(" - Mean unsigned error in all components: %14.7e\n", mue_f);
  printf(" - Mean unsigned error in the magnitude:  %14.7e\n", mue_f_mag);
  output = terminalFormat("Bond parameters represented in 32-bit floating point numbers, "
                          "coordinates represented in 64-bit signed integers, atomic "
                          "displacements computed as signed integers, all other calculations done "
                          "in 32-bit floating point arithmetic:", nullptr, nullptr, 0, 1, 1);
  printf("%s\n", output.c_str());
  printf(" - Mean unsigned error in all components: %14.7e\n", mue_flli);
  printf(" - Mean unsigned error in the magnitude:  %14.7e\n", mue_flli_mag);
  output = terminalFormat("Bond equilibrium represented as a 64-bit floating point number, "
                          "coordinates represented in 64-bit signed integers, atomic "
                          "displacements and computed as signed integers, bond stretch computed "
                          "in 64-bit arithmetic, all other calculations done in 32-bit floating "
                          "point arithmetic:", nullptr, nullptr, 0, 1, 1);
  printf("%s\n", output.c_str());
  printf(" - Mean unsigned error in all components: %14.7e\n", mue_f2lli);
  printf(" - Mean unsigned error in the magnitude:  %14.7e\n", mue_f2lli_mag);
  
  // Read topology
  const char osc = osSeparator();
  const std::string topology_home = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string trpcage_top = topology_home + osc + "trpcage.top";
  const std::string coordinate_home = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string trpcage_crd = coordinate_home + osc + "trpcage.inpcrd";
  AtomGraph trpcage_ag;
  PhaseSpace trpcage_ps;
  ScoreCard all_systems_sc(1);
  const int trpcage_idx = 0;
  if (getDrivePathType(trpcage_top) == DrivePathType::FILE) {
    trpcage_ag.buildFromPrmtop(trpcage_top, ExceptionResponse::SILENT);
  }
  else {
    rtErr("The Trp-cage topology " + trpcage_top + " was not found.  The benchmark cannot run.",
          "SplitForceAccumulation", "split_valence");
  }
  if (getDrivePathType(trpcage_crd) == DrivePathType::FILE) {
    trpcage_ps.buildFromFile(trpcage_crd, CoordinateFileKind::AMBER_INPCRD);
  }
  else {
    rtErr("The Trp-cage topology " + trpcage_crd + " was not found.  The benchmark cannot run.",
          "SplitForceAccumulation", "split_valence");
  }

  // Evaluate the forces due to bond and bond angle terms in double precision computations
  const TrajectoryKind tkind = TrajectoryKind::FORCES;
  trpcage_ps.initializeForces();
  const double trpcage_bond_e = evaluateBondTerms(trpcage_ag, &trpcage_ps, &all_systems_sc,
                                                  EvaluateForce::YES, trpcage_idx);
  const std::vector<double> trpcage_bond_frc = trpcage_ps.getInterlacedCoordinates(tkind);
  trpcage_ps.initializeForces();
  const double trpcage_angl_e = evaluateAngleTerms(trpcage_ag, &trpcage_ps, &all_systems_sc,
                                                   EvaluateForce::YES, trpcage_idx);
  const std::vector<double> trpcage_angl_frc = trpcage_ps.getInterlacedCoordinates(tkind);
  trpcage_ps.initializeForces();
}
