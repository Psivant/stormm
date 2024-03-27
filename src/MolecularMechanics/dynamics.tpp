// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace mm {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc, typename Tcalc2, typename Tcalc4>
void dynaStep(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, const Tcoord* xvel,
              const Tcoord* yvel, const Tcoord* zvel, Tcoord* xfrc, Tcoord* yfrc, Tcoord* zfrc,
              Tcoord* xalt, Tcoord* yalt, Tcoord* zalt, Tcoord* vxalt, Tcoord* vyalt,
              Tcoord* vzalt, Tcoord* fxalt, Tcoord* fyalt, Tcoord* fzalt, ScoreCard *sc,
              const ThermostatWriter<Tcalc> &tstw, const ValenceKit<Tcalc> &vk,
              const NonbondedKit<Tcalc> &nbk, const ImplicitSolventKit<Tcalc> &isk,
              const NeckGeneralizedBornKit<Tcalc> &neck_gbk, Tcoord* effective_gb_radii,
              Tcoord* psi, Tcoord* sumdeijda, const RestraintKit<Tcalc, Tcalc2, Tcalc4> &rar,
              const VirtualSiteKit<Tcalc> &vsk, const ChemicalDetailsKit &cdk,
              const ConstraintKit<Tcalc> &cnk, const StaticExclusionMaskReader &ser,
              const DynamicsControls &dyncon, const int system_index,
              const Tcalc gpos_scale_factor, const Tcalc vel_scale_factor,
              const Tcalc frc_scale_factor) {
  
  // Evaluate the force and energy for a system in vacuum with isolated boundary conditions
  evalRestrainedMMGB<Tcoord, Tcoord,
                     Tcalc, Tcalc2, Tcalc4>(xcrd, ycrd, zcrd, nullptr, nullptr, UnitCellType::NONE,
                                            xfrc, yfrc, zfrc, sc, vk, nbk, ser, isk, neck_gbk,
                                            effective_gb_radii, psi, sumdeijda, rar,
                                            EvaluateForce::YES, system_index, tstw.step);
  transmitVirtualSiteForces<Tcoord, Tcoord, Tcalc>(xcrd, ycrd, zcrd, xfrc, yfrc, zfrc, nullptr,
                                                   nullptr, UnitCellType::NONE, vsk);
  
  // Update the velocities by the first half step with the new forces.
  velocityVerletVelocityUpdate<Tcoord, Tcalc>(xvel, yvel, zvel, xfrc, yfrc, zfrc, cdk.natom,
                                              cdk.masses, vxalt, vyalt, vzalt, tstw,
                                              vel_scale_factor, frc_scale_factor);
  
  // Constrain velocities
  if (tstw.cnst_geom) {
    rattleVelocities<Tcoord, Tcalc>(vxalt, vyalt, vzalt, xcrd, ycrd, zcrd, cnk, tstw.dt,
                                    dyncon.getRattleTolerance(), dyncon.getRattleIterations(),
                                    dyncon.getCpuRattleMethod(), gpos_scale_factor,
                                    vel_scale_factor);
  }
  
  // Commit the energy, all components (energy computations are obligatory in CPU functions).  The
  // diagnostics from the initial state will always be stored.
  if (dyncon.getDiagnosticPrintFrequency() > 0 &&
      tstw.step % dyncon.getDiagnosticPrintFrequency() == 0) {
    evalKineticEnergy<Tcoord, Tcalc>(vxalt, vyalt, vzalt, sc, cdk, system_index,
                                     static_cast<Tcalc>(1.0) / vel_scale_factor);
    computeTemperature(sc, cdk, tstw.cnst_geom, system_index);
    sc->commit(StateVariable::ALL_STATES, system_index);
    sc->incrementSampleCount();
    sc->setLastTimeStep(tstw.step);
  }

  // Move particles, placing their new positions in the {x,y,z}alt arrays.
  velocityVerletCoordinateUpdate<Tcoord, Tcalc>(xcrd, ycrd, zcrd, xfrc, yfrc, zfrc, cdk.natom,
                                                cdk.masses, xalt, yalt, zalt, vxalt, vyalt, vzalt,
                                                tstw, gpos_scale_factor, vel_scale_factor,
                                                frc_scale_factor);

  // Apply positional constraints
  if (tstw.cnst_geom) {
    rattlePositions<Tcoord, Tcalc>(xalt, yalt, zalt, vxalt, vyalt, vzalt, xcrd, ycrd, zcrd, cnk,
                                   tstw.dt, dyncon.getRattleTolerance(),
                                   dyncon.getRattleIterations(), dyncon.getCpuRattleMethod(),
                                   gpos_scale_factor, vel_scale_factor);
  }

  // Replace virtual sites
  placeVirtualSites<Tcoord, Tcalc>(xalt, yalt, zalt, nullptr, nullptr, UnitCellType::NONE, vsk,
                                   gpos_scale_factor);

  // Zero forces in the alternate time point, in preparation for the next step.  Auxiliary arrays
  // involved in Generalized Born calculations (psi, effective_gb_radii, sumdeijda) will be
  // initialized in their respective CPU routines.
  const Tcoord zero = 0.0;
  for (int i = 0; i < cdk.natom; i++) {
    fxalt[i] = zero;
    fyalt[i] = zero;
    fzalt[i] = zero;
  }
}

} // namespace mm
} // namespace stormm
