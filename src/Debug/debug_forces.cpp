#include "copyright.h"
#include "Constants/behavior.h"
#include "MolecularMechanics/dynamics_intervention.h"
#include "MolecularMechanics/mm_evaluation.h"
#include "Potential/energy_enumerators.h"
#include "Potential/eval_synthesis.h"
#include "Trajectory/coordinate_copy.h"
#include "Trajectory/trajectory_enumerators.h"
#include "Topology/atomgraph_enumerators.h"
#include "debug.h"
#ifdef STORMM_USE_HPC
#include "hpc_debug.h"
#endif

namespace stormm {
namespace mm {

using constants::PrecisionModel;
using debug::checkAtomicForces;
using energy::evalSyNonbondedEnergy;
using energy::EvaluateForce;
using energy::NonbondedTask;
using topology::UnitCellType;
using trajectory::coordCopy;
using trajectory::TrajectoryKind;

//-------------------------------------------------------------------------------------------------
void execForceDebug(const int step, const int index, const HybridTargetLevel tier) {
  const PhaseSpaceSynthesis *poly_ps = dyna_tk.getPhaseSpaceSynthesisPointer();
  PhaseSpaceSynthesis *workspc = dyna_tk.getWorkspacePointer();
  AtomGraphSynthesis *poly_ag = dyna_tk.getTopologySynthesisPointer();
  const PsSynthesisReader poly_psr = poly_ps->data(tier);
  PsSynthesisWriter workspc_w = workspc->data();
  const GpuDetails gpu = dyna_tk.getGpuDetails();
  
  // If the coordinate synthesis is not prepared to accept a download of coordinates, that will
  // become clear at this step.  The coordinates are downloaded one ysstem at a time, so that only
  // positions for the active point in the WHITE / BLACK coordinate cycle can be transferred
  // (rather than also transferring velocities and forces).
  const int nsys = poly_ps->getSystemCount();
  for (int i = 0; i < nsys; i++) {
    coordCopy(&workspc_w, workspc_w.atom_starts[i], i, poly_psr, workspc_w.atom_starts[i], i,
              workspc_w.atom_counts[i], HybridTargetLevel::HOST, tier, gpu);
  }
  ScoreCard sc(poly_ps->getSystemCount(), 1, 36);

  // Initialize forces on the CPU host, then compute them for each system's configuration.
  workspc->initializeForces();
  switch (poly_ps->getUnitCellType()) {
  case UnitCellType::NONE:
    {
      // Compute valence forces for all systems
      evalValeRestMM(workspc, &sc, poly_ag, step, EvaluateForce::YES, PrecisionModel::DOUBLE);

      // Compute non-bonded interactions for all systems
      const StaticExclusionMaskSynthesis *syse = dyna_tk.getStaticExclusionMaskPointer();
      if (dyna_tk.hasImplicitSolventWorkspace()) {
        ImplicitSolventWorkspace *isw_ptr = dyna_tk.getImplicitSolventWorkspacePointer();
        isw_ptr->initialize();
        evalSyNonbondedEnergy(*poly_ag, *syse, workspc, isw_ptr, &sc, NonbondedTask::GB_RADII,
                              PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateForce::YES);
        evalSyNonbondedEnergy(*poly_ag, *syse, workspc, isw_ptr, &sc,
                              NonbondedTask::PARTICLE_PARTICLE, PrecisionModel::DOUBLE,
                              EvaluateForce::YES, EvaluateForce::YES);
        evalSyNonbondedEnergy(*poly_ag, *syse, workspc, isw_ptr, &sc,
                              NonbondedTask::GB_RADII_DERIVATIVES, PrecisionModel::DOUBLE,
                              EvaluateForce::YES, EvaluateForce::YES);      
      }
      else {
        ImplicitSolventWorkspace isw(poly_ag->getSystemAtomOffsets(),
                                     poly_ag->getSystemAtomCounts(), PrecisionModel::DOUBLE);
        isw.initialize();
        evalSyNonbondedEnergy(*poly_ag, *syse, workspc, &isw, &sc, NonbondedTask::GB_RADII,
                              PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateForce::YES);
        evalSyNonbondedEnergy(*poly_ag, *syse, workspc, &isw, &sc,
                              NonbondedTask::PARTICLE_PARTICLE, PrecisionModel::DOUBLE,
                              EvaluateForce::YES, EvaluateForce::YES);
        evalSyNonbondedEnergy(*poly_ag, *syse, workspc, &isw, &sc,
                              NonbondedTask::GB_RADII_DERIVATIVES, PrecisionModel::DOUBLE,
                              EvaluateForce::YES, EvaluateForce::YES);      
      }
    }  
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    break;
  }
  
  // Compare forces just computed on the CPU host in the workspace against those in the original
  // calculations.
  for (int i = 0; i < workspc_w.system_count; i++) {
    const CoordinateFrame orig_cf = poly_ps->exportCoordinates(i, HybridFormat::HOST_ONLY,
                                                               TrajectoryKind::FORCES, tier);
    const CoordinateFrame chek_cf = workspc->exportCoordinates(i, HybridFormat::HOST_ONLY,
                                                               TrajectoryKind::FORCES,
                                                               HybridTargetLevel::HOST);
    const CoordinateFrameReader orig_cfr = orig_cf.data();
    const CoordinateFrameReader chek_cfr = orig_cf.data();
    for (int j = 0; j < orig_cfr.natom; j++) {

      // CHECK
      if (fabs(orig_cfr.xcrd[j] - chek_cfr.xcrd[j]) > 5.0e-3 ||
          fabs(orig_cfr.ycrd[j] - chek_cfr.ycrd[j]) > 5.0e-3 ||
          fabs(orig_cfr.zcrd[j] - chek_cfr.zcrd[j]) > 5.0e-3) {
        printf("Atom %4d : %6d  Step %6d    %9.4lf %9.4lf %9.4lf    %9.4lf %9.4lf %9.4lf\n", i, j,
               step, orig_cfr.xcrd[j], orig_cfr.ycrd[j], orig_cfr.zcrd[j], chek_cfr.xcrd[j],
               chek_cfr.ycrd[j], chek_cfr.zcrd[j]);
        exit(1);
      }
      // END CHECK
    }
  }
}

//-------------------------------------------------------------------------------------------------
void execNLForceDebug(const int step, const int index, const HybridTargetLevel tier) {

}

//-------------------------------------------------------------------------------------------------
void evalForceAnomalies(const int step, const int index, const HybridTargetLevel tier) {
  PsSynthesisReader poly_psr = dyna_tk.getReadOnlyCoordinateData(tier);
  Watcher *bug = dyna_tk.getAnomalyReportingPointer();
  WatcherWriter bugw = dyna_tk.getAnomalyReportingData(tier);
  const int init_nfrc_events = bug->getLargeForceCount(tier);
  switch (tier) {
  case HybridTargetLevel::HOST:
    checkAtomicForces(&bugw, poly_psr, step, IntegrationStage::CALC_FORCES);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    checkAtomicForces(&bugw, poly_psr, step, IntegrationStage::CALC_FORCES,
                      dyna_tk.getGpuDetails());
    break;
#endif
  }
  const int finl_nfrc_events = bug->getLargeForceCount(tier);
  if (finl_nfrc_events != init_nfrc_events) {
    execForceDebug(step, index, tier);
  }
}
  
} // namespace mm
} // namespace stormm
