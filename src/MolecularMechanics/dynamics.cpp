#include "copyright.h"
#include "FileManagement/file_enumerators.h"
#include "Trajectory/trajectory_enumerators.h"
#include "Trajectory/trim.h"
#include "dynamics.h"

namespace stormm {
namespace mm {

using trajectory::CoordinateFileKind;
using trajectory::TrajectoryKind;
using diskutil::PrintSituation;
  
//-------------------------------------------------------------------------------------------------
void dynaStep(PhaseSpaceWriter *psw, ScoreCard *sc, const ThermostatWriter<double> &tstw,
              const ValenceKit<double> &vk, const NonbondedKit<double> &nbk,
              const ImplicitSolventKit<double> &isk,
              const NeckGeneralizedBornKit<double> &neck_gbk, double* effective_gb_radii,
              double* psi, double* sumdeijda, const RestraintKit<double, double2, double4> &rar,
              const VirtualSiteKit<double> &vsk, const ChemicalDetailsKit &cdk,
              const ConstraintKit<double> &cnk, const StaticExclusionMaskReader &ser,
              const DynamicsControls &dyncon, const int system_index) {
  dynaStep<double, double,
           double2, double4>(psw->xcrd, psw->ycrd, psw->zcrd, psw->xvel, psw->yvel, psw->zvel,
                             psw->xfrc, psw->yfrc, psw->zfrc, psw->xalt, psw->yalt, psw->zalt,
                             psw->vxalt, psw->vyalt, psw->vzalt, psw->fxalt, psw->fyalt,
                             psw->fzalt, sc, tstw, vk, nbk, isk, neck_gbk, effective_gb_radii, psi,
                             sumdeijda, rar, vsk, cdk, cnk, ser, dyncon, system_index);
}
  
//-------------------------------------------------------------------------------------------------
void dynamics(PhaseSpace *ps, Thermostat *heat_bath, ScoreCard *sc, const AtomGraph *ag,
              const NeckGeneralizedBornTable *neck_gbtab, const StaticExclusionMask *se,
              const RestraintApparatus *ra, const DynamicsControls &dyncon,
              const int system_index, const std::string &trajectory_file_name,
              const std::string &restart_file_name) {

  // Produce abstracts at each point in the coordinate object's time cycle.
  PhaseSpaceWriter psw = ps->data();
  PhaseSpaceWriter psw_alt = ps->data(getNextCyclePosition(ps->getCyclePosition()));

  // Produce critical abstracts for other objects as needed.
  ThermostatWriter tstw = heat_bath->dpData();
  const ValenceKit<double> vk = ag->getDoublePrecisionValenceKit();
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  const ImplicitSolventKit<double> isk = ag->getDoublePrecisionImplicitSolventKit();
  const VirtualSiteKit<double> vsk = ag->getDoublePrecisionVirtualSiteKit();
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  const ConstraintKit<double> cnk = ag->getDoublePrecisionConstraintKit();
  const NeckGeneralizedBornKit<double> neck_gbk = neck_gbtab->dpData();
  const RestraintKit<double, double2, double4> rar = ra->dpData();
  const StaticExclusionMaskReader ser = se->data();

  // Store critical output constants
  const int traj_freq = dyncon.getTrajectoryPrintFrequency();
  const int cmpg_freq = dyncon.getCenterOfMassMotionPurgeFrequency();
  
  // Allocate arrays to perform GB calculations.
  std::vector<double> effective_gb_radii(psw.natom), psi(psw.natom), sumdeijda(psw.natom);
  
  // Loop for the requested number of cycles.  One calculation alternates between two abstracts
  // taking the coordinates at alternating points in the time cycle.
  for (int step = 0; step < dyncon.getStepCount(); step++) {

    // If the thermostat's random number cache has, by this time, been used up, refresh it.
    if (step % tstw.depth == 0 && step > 0) {
      heat_bath->refresh(0, ps->getAtomCount());
    }

    // If requested, purge motion of the center of mass for each system.
    if (cmpg_freq > 0 && step > 0 && step % cmpg_freq == 0) {
      removeMomentum(ps, ag);
    }
    
    // Compute forces, update velocities, and move particles in the velocity-Verlet scheme
    if (step & 0x1) {
      dynaStep(&psw_alt, sc, tstw, vk, nbk, isk, neck_gbk, effective_gb_radii.data(), psi.data(),
               sumdeijda.data(), rar, vsk, cdk, cnk, ser, dyncon, system_index);
    }
    else {
      dynaStep(&psw, sc, tstw, vk, nbk, isk, neck_gbk, effective_gb_radii.data(), psi.data(),
               sumdeijda.data(), rar, vsk, cdk, cnk, ser, dyncon, system_index);
    }

    // If trajectory frames are required, produce one now.  State the time in units of ps to be
    // consistent with Amber.
    if (traj_freq > 0 && (step + 1) % traj_freq == 0) {
      const double current_time = static_cast<double>(step) * tstw.dt;
      ps->exportToFile(trajectory_file_name, current_time, TrajectoryKind::POSITIONS,
                       CoordinateFileKind::AMBER_CRD, PrintSituation::APPEND);
    }
    
    // Advance the coordinate object's time cycle.  This is done for consistency and future
    // reference: the correct abstract will be chosen based on the modulo operation above.
    ps->updateCyclePosition();
    
    // Update the time step.  The abstract can be updated without completely regenerating it.
    tstw.step += 1;
    heat_bath->incrementStep();
  }

  // Total up the potential and total energies
  sc->computePotentialEnergy();
  sc->computeTotalEnergy();
}

//-------------------------------------------------------------------------------------------------
void dynamics(PhaseSpace *ps, Thermostat *heat_bath, ScoreCard *sc, const AtomGraph &ag,
              const NeckGeneralizedBornTable &neck_gbtab, const StaticExclusionMask &se,
              const RestraintApparatus &ra, const DynamicsControls &dyncon,
              const int system_index, const std::string &trajectory_file_name,
              const std::string &restart_file_name) {
  dynamics(ps, heat_bath, sc, ag.getSelfPointer(), neck_gbtab.getSelfPointer(),
           se.getSelfPointer(), ra.getSelfPointer(), dyncon, system_index, trajectory_file_name,
           restart_file_name);
}

} // namespace mm
} // namespace stormm
