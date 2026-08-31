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

  // Find the mass array that is amenable to the template.  This is not an ideal thing to do, but
  // it will cast the float* array to Tcalc* if Tcalc is float and the double* array to Tcalc* if
  // Tcalc is double.  This requirement is the legacy of ChemicalDetailsKit originally not
  // containing the array of masses and therefore having no need to be templated in its own right.
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const Tcalc* mass_ptr = (tcalc_is_double) ? reinterpret_cast<const Tcalc*>(cdk.masses) :
                                              reinterpret_cast<const Tcalc*>(cdk.sp_masses); 
 
  // Update the velocities by the first half step with the new forces.
  velocityVerletVelocityUpdate<Tcoord, Tcalc>(xvel, yvel, zvel, xfrc, yfrc, zfrc, cdk.natom,
                                              mass_ptr, vxalt, vyalt, vzalt, tstw, nullptr,
                                              nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                                              nullptr, nullptr, 0, vel_scale_factor,
                                              frc_scale_factor);
  
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
                                                mass_ptr, xalt, yalt, zalt, vxalt, vyalt, vzalt,
                                                tstw, nullptr, nullptr, nullptr, nullptr, nullptr,
                                                nullptr, nullptr, nullptr, nullptr, nullptr,
                                                nullptr, nullptr, 0, gpos_scale_factor,
                                                vel_scale_factor, frc_scale_factor);

  // Apply positional constraints
  if (tstw.cnst_geom) {
    shakePositions<Tcoord, Tcalc>(xalt, yalt, zalt, vxalt, vyalt, vzalt, xcrd, ycrd, zcrd, cnk,
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
  
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcoord4, typename Tval_calc,
          typename Tval_calc2, typename Tval_calc4, typename Tnb_calc, typename Tnb_calc2>
void dynaStep(PsSynthesisWriter *poly_psw, const PsSynthesisBorders &pssb,
              CellGridWriter<void, void, void, void> *cgw_v, PMIGridAccumulator *pm_acc,
              PMIGridWriter *pm_wrt, const PMIGridReader &pm_rdr,
              ConvolutionWriter<Tnb_calc, Tnb_calc2> *cvolw, ScoreCard *sc,
              ThermostatWriter<Tval_calc> *tstw, const SyValenceKit<Tval_calc> &poly_vk,
              const SyNonbondedKit<Tnb_calc, Tnb_calc2> &poly_nbk,
              const SyRestraintKit<Tval_calc, Tval_calc2, Tval_calc4> &poly_rk,
              const SyAtomUpdateKit<Tval_calc, Tval_calc2, Tval_calc4> &poly_auk,
              const LocalExclusionMaskReader &lemr, const Tnb_calc cutoff,
              const Tnb_calc qqew_coeff, const VdwSumMethod vdw_sum, const int ntpr,
              const int ntwx) {
  
  // Create read-only forms of the coordinate synthesis and thermostat abstracts
  const PsSynthesisReader poly_psr(poly_psw);
  const ThermostatReader<Tval_calc> tstr(tstw);

  // Produce a mutable energy tracking abstract
  ScoreCardWriter scw = sc->data();
  
  // Compute the non-bonded particle-particle interactions
  CellGridWriter<Tcoord, Tacc,
                 Tnb_calc, Tcoord4> cgw = restoreType<Tcoord, Tacc, Tnb_calc, Tcoord4>(cgw_v);
  const CellGridReader<Tcoord, Tacc, Tnb_calc, Tcoord4> cgr(cgw);
  evaluateParticleParticleEnergy<Tcoord, Tacc,
                                 Tnb_calc, Tnb_calc2, Tcoord4>(cgw_v, &scw, poly_psr, poly_nbk,
                                                               lemr, cutoff, qqew_coeff, vdw_sum,
                                                               EvaluateForce::YES,
                                                               NonbondedTheme::ALL);

  // Compute the non-bonded particle-mesh interactions
  evaluateParticleMeshEnergy<Tcoord, Tacc,
                             Tnb_calc, Tnb_calc2, Tcoord4>(&cgw, cvolw, pm_acc, pm_wrt, pm_rdr,
                                                           pssb, poly_nbk, &scw);

  // Contribute the non-bonded forces to the synthesis accumulators
  contributeCellGridForces<Tcoord, Tacc, Tnb_calc, Tcoord4>(poly_psw, cgr);

  // Compute valence interactions
  evalValeRestMM<Tval_calc, Tval_calc2, Tval_calc4>(poly_psw, sc, poly_vk, poly_rk, poly_auk,
                                                    EvaluateForce::YES, VwuTask::ALL_TASKS,
                                                    tstw->step);
  
  // Transmit virtual site forces to frame atoms, if applicable
  transmitVirtualSiteForces<Tval_calc, Tval_calc2, Tval_calc4>(poly_psw, poly_vk, poly_auk);  
  
  // Velocity Verlet update I
  velocityVerletVelocityUpdate<Tval_calc, Tval_calc2, Tval_calc4>(poly_psw, poly_auk, tstr);

  // Apply velocity constraints
  if (tstr.cnst_geom) {
    rattleVelocities<Tval_calc, Tval_calc2, Tval_calc4>(poly_psw, poly_vk, poly_auk, tstr.dt,
                                                        tstr.rattle_tol, tstr.rattle_iter);
    settleVelocities<Tval_calc, Tval_calc2, Tval_calc4>(poly_psw, poly_vk, poly_auk);
  }

  // Calculate kinetic energy and temperature in step with evaluations of the potential energy
  if (ntpr > 0 && tstw->step % ntpr == 0) {
    evalKineticEnergy<Tval_calc, Tval_calc2, Tval_calc4, Tval_calc>(poly_psr, sc, poly_auk);
    computeTemperature<Tval_calc, Tval_calc2, Tval_calc4, Tval_calc>(poly_psr, sc, poly_auk, tstr,
                                                                     true);
    sc->commit(StateVariable::ALL_STATES);
    sc->incrementSampleCount();
    sc->setLastTimeStep(tstw->step);
  }
  
  // Velocity Verlet update II
  velocityVerletCoordinateUpdate<Tval_calc, Tval_calc2, Tval_calc4>(poly_psw, poly_auk, tstr);
  
  // Apply position constraints
  if (tstr.cnst_geom) {
    shakePositions<Tval_calc>(poly_psw, poly_vk, poly_auk, tstr.dt, tstr.rattle_tol,
                              tstr.rattle_iter);
    settlePositions<Tval_calc, Tval_calc2, Tval_calc4>(poly_psw, poly_vk, poly_auk, tstr.dt);
  }

  // Replace virtual sites
  placeVirtualSites<Tval_calc, Tval_calc2, Tval_calc4>(poly_psw, poly_vk, poly_auk);

  // Interventions intended to occur after geometry constraints will also take place on the
  // repositioned virtual particles.
  
  // Migrate particles within the cell grid
  migrate(&cgw, poly_psr);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void dynamics(PhaseSpaceSynthesis *poly_ps, CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> *cg,
              PMIGrid *pmig, ConvolutionManager *cvol, ScoreCard *sc, Thermostat *heat_bath,
              const AtomGraphSynthesis &poly_ag, const LocalExclusionMask &lem,
              const DynamicsControls &dyncon, const PrecisionControls &preccon,
              const PPPMControls &pmecon) {

  // Check that the CellGrid calculations match the stated non-bonded calculation type.
  bool calc_mismatch = false;
  switch (preccon.getNonbondedMethod()) {
  case PrecisionModel::DOUBLE:
    if (std::type_index(typeid(Tcalc)).hash_code() == float_type_index) {
      calc_mismatch = true;
    }
    break;
  case PrecisionModel::SINGLE:
    if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
      calc_mismatch = true;
    }
    break;
  }
  if (calc_mismatch) {
    rtErr("Mismatched types were submitted for the CellGrid calculations (" +
          getStormmScalarTypeName<Tcalc>() + ") and the general non-bonded methods (" +
          getEnumerationName(preccon.getNonbondedMethod()) + ").", "dynamics");
  }
  
  // Produce abstracts of the relevant objects, at each point in the coordinate time cycle
  const CoordinateCycle poly_ps_next_stage = getNextCyclePosition(poly_ps->getCyclePosition());
  PsSynthesisWriter poly_psw = poly_ps->data();
  PsSynthesisWriter poly_psw_alt = poly_ps->data(poly_ps_next_stage);
  const PsSynthesisBorders pssb = poly_ps->borders();
  const PsSynthesisBorders pssb_alt = poly_ps->borders(poly_ps_next_stage);
  ScoreCardWriter scw = sc->data();
  const LocalExclusionMaskReader lemr = lem.data();
  const CoordinateCycle cg_next_stage = getNextCyclePosition(cg->getCyclePosition());
  CellGridWriter<void, void, void, void> cgv = cg->templateFreeData();
  CellGridWriter<void, void, void, void> cgv_alt = cg->templateFreeData(cg_next_stage);
  ConvolutionWriter<double, double2> d_cvolw = cvol->dpData();
  ConvolutionWriter<float, float2> f_cvolw = cvol->spData();
  PMIGridAccumulator pm_acc = pmig->fpData(HybridTargetLevel::HOST, ExceptionResponse::SILENT);
  PMIGridWriter pm_wrt = pmig->data();
  const PMIGridReader pm_rdr(pm_wrt);
  const SyNonbondedKit<double, double2> dpoly_nbk = poly_ag.getDoublePrecisionNonbondedKit();
  const SyNonbondedKit<float, float2> fpoly_nbk = poly_ag.getSinglePrecisionNonbondedKit();
  MotionSweeper mos(poly_ps);
  ThermostatWriter<double> d_tstw = heat_bath->dpData();
  ThermostatWriter<float> f_tstw = heat_bath->spData();
  const SyValenceKit<double> dpoly_vk = poly_ag.getDoublePrecisionValenceKit();
  const SyValenceKit<float> fpoly_vk = poly_ag.getSinglePrecisionValenceKit();
  const SyAtomUpdateKit<double,
                        double2,
                        double4_16a> dpoly_auk = poly_ag.getDoublePrecisionAtomUpdateKit();
  const SyAtomUpdateKit<float,
                        float2, float4> fpoly_auk = poly_ag.getSinglePrecisionAtomUpdateKit();
  const SyRestraintKit<double,
                       double2, double4_16a> dpoly_rk = poly_ag.getDoublePrecisionRestraintKit();
  const SyRestraintKit<float,
                       float2, float4> fpoly_rk = poly_ag.getSinglePrecisionRestraintKit();
  const int traj_freq = dyncon.getTrajectoryPrintFrequency();
  const int diag_freq = dyncon.getDiagnosticPrintFrequency();
  const int cmpg_freq = dyncon.getCenterOfMassMotionPurgeFrequency();

  // Cutoffs and critical non-bonded constants.  Selection of the Ewald coefficient would be done
  // for the setup of objects used in the GPU workflow, but on the CPU the Ewald coefficient can be
  // extracted 
  const double cutoff = cg->getCutoff();
  const double dsum_tol = pmecon.getDirectSumTolerance();
  const double ew_coeff = ewaldCoefficient(cutoff, dsum_tol);
  const VdwSumMethod vdw_sum = dyncon.getVdwSummation();
  
  // Produce critical abstracts for other objects as needed
  for (int step = 0; step < dyncon.getStepCount(); step++) {
    
    // If the thermostat's random number cache has, by this time, been used up, refresh it.  While
    // one one of the thermostat abstracts is used for the intensive dynamics calculations, both
    // contain valid information about the random number cache depth and the overall padded atom
    // count for the synthesis.
    if (step > 0 && d_tstw.depth > 0 && step % d_tstw.depth == 0) {
      heat_bath->refresh(0, d_tstw.padded_natom);
    }

    // If requested, purge motion of the center of mass for each system.
    if (cmpg_freq > 0 && step > 0 && step % cmpg_freq == 0) {
      removeMomentum(poly_ps, poly_ag, &mos);
    }

    // Initialize forces and energy accumulators.
    poly_ps->initializeForces();
    cg->initializeForces();
    pmig->initialize();
    sc->initialize();
        
    // Perform the step.
    if (step & 0x1) {
      switch (preccon.getValenceMethod()) {
      case PrecisionModel::DOUBLE:
        switch (preccon.getNonbondedMethod()) {
        case PrecisionModel::DOUBLE:
          dynaStep<Tcoord, Tacc,
                   Tcoord4, double,
                   double2, double4_16a,
                   double, double2>(&poly_psw_alt, pssb_alt, &cgv_alt, &pm_acc, &pm_wrt, pm_rdr,
                                    &d_cvolw, sc, &d_tstw, dpoly_vk, dpoly_nbk, dpoly_rk,
                                    dpoly_auk, lemr, cutoff, ew_coeff, vdw_sum, diag_freq,
                                    traj_freq);
          break;
        case PrecisionModel::SINGLE:
          dynaStep<Tcoord, Tacc,
                   Tcoord4, double,
                   double2, double4_16a,
                   float, float2>(&poly_psw_alt, pssb_alt, &cgv_alt, &pm_acc, &pm_wrt, pm_rdr,
                                  &f_cvolw, sc, &d_tstw, dpoly_vk, fpoly_nbk, dpoly_rk, dpoly_auk,
                                  lemr, cutoff, ew_coeff, vdw_sum, diag_freq, traj_freq);
          break;
        }
        break;
      case PrecisionModel::SINGLE:
        switch (preccon.getNonbondedMethod()) {
        case PrecisionModel::DOUBLE:
          dynaStep<Tcoord, Tacc,
                   Tcoord4, float,
                   float2, float4,
                   double, double2>(&poly_psw_alt, pssb_alt, &cgv_alt, &pm_acc, &pm_wrt, pm_rdr,
                                    &d_cvolw, sc, &f_tstw, fpoly_vk, dpoly_nbk, fpoly_rk,
                                    fpoly_auk, lemr, cutoff, ew_coeff, vdw_sum, diag_freq,
                                    traj_freq);
          break;
        case PrecisionModel::SINGLE:
          dynaStep<Tcoord, Tacc,
                   Tcoord4, float,
                   float2, float4,
                   float, float2>(&poly_psw_alt, pssb_alt, &cgv_alt, &pm_acc, &pm_wrt, pm_rdr,
                                  &f_cvolw, sc, &f_tstw, fpoly_vk, fpoly_nbk, fpoly_rk, fpoly_auk,
                                  lemr, cutoff, ew_coeff, vdw_sum, diag_freq, traj_freq);
          break;
        }
        break;
      }
    }
    else {
      switch (preccon.getValenceMethod()) {
      case PrecisionModel::DOUBLE:
        switch (preccon.getNonbondedMethod()) {
        case PrecisionModel::DOUBLE:
          dynaStep<Tcoord, Tacc,
                   Tcoord4, double,
                   double2, double4_16a,
                   double, double2>(&poly_psw, pssb, &cgv, &pm_acc, &pm_wrt, pm_rdr, &d_cvolw, sc,
                                    &d_tstw, dpoly_vk, dpoly_nbk, dpoly_rk, dpoly_auk, lemr,
                                    cutoff, ew_coeff, vdw_sum, diag_freq, traj_freq);
          break;
        case PrecisionModel::SINGLE:
          dynaStep<Tcoord, Tacc,
                   Tcoord4, double,
                   double2, double4_16a,
                   float, float2>(&poly_psw, pssb, &cgv, &pm_acc, &pm_wrt, pm_rdr, &f_cvolw, sc,
                                  &d_tstw, dpoly_vk, fpoly_nbk, dpoly_rk, dpoly_auk, lemr, cutoff,
                                  ew_coeff, vdw_sum, diag_freq, traj_freq);
          break;
        }
        break;
      case PrecisionModel::SINGLE:
        switch (preccon.getNonbondedMethod()) {
        case PrecisionModel::DOUBLE:
          dynaStep<Tcoord, Tacc,
                   Tcoord4, float,
                   float2, float4,
                   double, double2>(&poly_psw, pssb, &cgv, &pm_acc, &pm_wrt, pm_rdr, &d_cvolw, sc,
                                    &f_tstw, fpoly_vk, dpoly_nbk, fpoly_rk, fpoly_auk, lemr,
                                    cutoff, ew_coeff, vdw_sum, diag_freq, traj_freq);
          break;
        case PrecisionModel::SINGLE:
          dynaStep<Tcoord, Tacc,
                   Tcoord4, float,
                   float2, float4,
                   float, float2>(&poly_psw, pssb, &cgv, &pm_acc, &pm_wrt, pm_rdr, &f_cvolw, sc,
                                  &f_tstw, fpoly_vk, fpoly_nbk, fpoly_rk, fpoly_auk, lemr, cutoff,
                                  ew_coeff, vdw_sum, diag_freq, traj_freq);
          break;
        }
        break;
      }
    }

    // Update the cycle positions and time step.
    poly_ps->updateCyclePosition();
    cg->updateCyclePosition();
    d_tstw.step += 1;
    f_tstw.step += 1;
    heat_bath->incrementStep();
  }

  // Total up the potential and total energies
  sc->computePotentialEnergy();
  sc->computeTotalEnergy();
}

} // namespace mm
} // namespace stormm
