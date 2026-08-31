// -*-c++-*-
#ifndef STORMM_DYNAMICS_LIBRARY_H
#define STORMM_DYNAMICS_LIBRARY_H

#include "copyright.h"
#include "Constants/generalized_born.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "MolecularMechanics/mm_evaluation.h"
#include "Namelists/nml_dynamics.h"
#include "Namelists/nml_pppm.h"
#include "Namelists/nml_precision.h"
#include "Potential/cellgrid.h"
#include "Potential/convolution_manager.h"
#include "Potential/energy_enumerators.h"
#include "Potential/map_density.h"
#include "Potential/pme_potential.h"
#include "Potential/pme_util.h"
#include "Potential/pmigrid.h"
#include "Potential/scorecard.h"
#include "Potential/local_exclusionmask.h"
#include "Potential/static_exclusionmask.h"
#include "Restraints/restraint_apparatus.h"
#include "Structure/hub_and_spoke.h"
#include "Structure/settle.h"
#include "Structure/virtual_site_handling.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/integration.h"
#include "Trajectory/motion_sweeper.h"
#include "Trajectory/phasespace.h"
#include "Trajectory/thermostat.h"
#include "Trajectory/trajectory_enumerators.h"
#include "Trajectory/trim.h"
#include "kinetic.h"

namespace stormm {
namespace mm {

using data_types::getStormmScalarTypeName;
using energy::CellGrid;
using energy::CellGridReader;
using energy::CellGridWriter;
using energy::ConvolutionManager;
using energy::ConvolutionWriter;
using energy::contributeCellGridForces;
using energy::evaluateParticleMeshEnergy;
using energy::evaluateParticleParticleEnergy;
using energy::ewaldCoefficient;
using energy::LocalExclusionMask;
using energy::LocalExclusionMaskReader;
using energy::NonbondedTheme;
using energy::PMIGrid;
using energy::PMIGridAccumulator;
using energy::PMIGridReader;
using energy::PMIGridWriter;
using energy::restoreType;
using energy::ScoreCard;
using energy::ScoreCardWriter;
using energy::StaticExclusionMask;
using energy::StaticExclusionMaskReader;
using energy::VdwSumMethod;
using namelist::DynamicsControls;
using namelist::PrecisionControls;
using namelist::PPPMControls;
using restraints::RestraintApparatus;
using restraints::RestraintKit;
using structure::placeVirtualSites;
using structure::rattleVelocities;
using structure::settlePositions;
using structure::settleVelocities;
using structure::shakePositions;
using structure::transmitVirtualSiteForces;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisReader;
using synthesis::PsSynthesisWriter;
using synthesis::SyAtomUpdateKit;
using synthesis::SyNonbondedKit;
using synthesis::SyRestraintKit;
using synthesis::SyValenceKit;
using topology::AtomGraph;
using topology::ConstraintKit;
using topology::ImplicitSolventKit;
using topology::NonbondedKit;
using topology::ValenceKit;
using topology::VirtualSiteKit;
using trajectory::IntegrationStage;
using trajectory::MotionSweeper;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;
using trajectory::removeMomentum;
using trajectory::Thermostat;
using trajectory::ThermostatReader;
using trajectory::ThermostatWriter;
using trajectory::velocityVerletVelocityUpdate;
using trajectory::velocityVerletCoordinateUpdate;
using namespace generalized_born_defaults;
  
/// \brief Run the canonical MD dynamics step for systems in implicit or explicit solvent.
///
/// \param xcrd                Current positions of particles along the Cartesian X axis
/// \param ycrd                Current positions of particles along the Cartesian Y axis
/// \param zcrd                Current positions of particles along the Cartesian Z axis
/// \param xvel                Current velocities of particles in the Cartesian X direction
/// \param yvel                Current velocities of particles in the Cartesian Y direction
/// \param zvel                Current velocities of particles in the Cartesian Z direction
/// \param xfrc                Forces acting on all particles in the Cartesian X direction.  These
///                            forces will be accumulated in one step and then used to move
///                            particles.
/// \param yfrc                Forces acting on all particles in the Cartesian Y direction
/// \param zfrc                Forces acting on all particles in the Cartesian Z direction
/// \param xalt                Developing positions of particles in the Cartesian X direction
/// \param yalt                Developing positions of particles in the Cartesian Y direction
/// \param zalt                Developing positions of particles in the Cartesian Z direction
/// \param vxalt               Developing velocities of particles in the Cartesian X direction
/// \param vyalt               Developing velocities of particles in the Cartesian Y direction
/// \param vzalt               Developing velocities of particles in the Cartesian Z direction
/// \param fxalt               Alternate array for forces acting on particles in the Cartesian X
///                            direction (these will be initialized on each step in preparation for
///                            accumulations in the subsequent step)
/// \param fyalt               Alternate forces acting on particles in the Cartesian Y direction
/// \param fzalt               Alternate forces acting on particles in the Cartesian Z direction
/// \param psw                 Abstract of the coordinates, containing positions, velocities, and
///                            and forces for particles in one system, in both current and future
///                            ("alt") storage arrays
/// \param tstr                Abstract of the thermostat object (see description below)
/// \param tst                 Thermostat (integrator) containing the time step and various
///                            random number generators for maintaining temperature control
/// \param vk                  Valence interaction parameters for all bonded terms
/// \param nbk                 Non-bonded parameters for all atoms in the simulation
/// \param isk                 Implicit solvent parameters and thread-block specific resources
/// \param neck_gbk            "Neck" Generalized Born parameter tables
/// \param effective_gb_radii  Array to accumulate the effective GB radii for all atoms
/// \param psi                 Array to hold the GB Psi values for all atoms
/// \param sumdeijda           Array to accumulate the GB derivative factors for all atoms
/// \param rar                 Restraint term parameters for the system
/// \param vsk                 Parameters and frame atom indices for virtual sites in the system
/// \param cnk                 Information on constraints in the system
/// \param ser                 Static exclusion mask abstract containing information on which
///                            non-bonded interactions are omitted
/// \param dyncon              User (or developer-mocked) input from a &dynamics namelist block
/// \param nrg_scale_bits      The number of bits after the point of any energy component
///                            representation
/// \param gpos_scale_factor   Scaling factor on the global positions fixed-point representation
/// \param vel_scale_factor    Scaling factor on the particle velocities fixed-point representation
/// \param frc_scale_factor    Scaling factor on fixed-point force accumulation
/// \{
template <typename Tcoord, typename Tcalc, typename Tcalc2, typename Tcalc4>
void dynaStep(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, const Tcoord* xvel,
              const Tcoord* yvel, const Tcoord* zvel, Tcoord* xfrc, Tcoord* yfrc, Tcoord* zfrc,
              Tcoord* xalt, Tcoord* yalt, Tcoord* zalt, Tcoord* vxalt, Tcoord* vyalt,
              Tcoord* vzalt, Tcoord* fxalt, Tcoord* fyalt, Tcoord* fzalt, ScoreCard *sc,
              const ThermostatWriter<Tcalc> &tstr, const ValenceKit<Tcalc> &vk,
              const NonbondedKit<Tcalc> &nbk, const ImplicitSolventKit<Tcalc> &isk,
              const NeckGeneralizedBornKit<Tcalc> &neck_gbk, Tcoord* effective_gb_radii,
              Tcoord* psi, Tcoord* sumdeijda, const RestraintKit<Tcalc, Tcalc2, Tcalc4> &rar,
              const VirtualSiteKit<Tcalc> &vsk, const ChemicalDetailsKit &cdk,
              const ConstraintKit<Tcalc> &cnk, const StaticExclusionMaskReader &ser,
              const DynamicsControls &dyncon, int system_index = 0, Tcalc gpos_scale_factor = 1.0,
              Tcalc vel_scale_factor = 1.0, Tcalc frc_scale_factor = 1.0);

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
              const LocalExclusionMaskReader &lemr, Tnb_calc cutoff, Tnb_calc qqew_coeff,
              VdwSumMethod vdw_sum, int ntpr, int ntwx);

void dynaStep(PhaseSpaceWriter *psw, ScoreCard *sc, const ThermostatWriter<double> &tstr,
              const ValenceKit<double> &vk, const NonbondedKit<double> &nbk,
              const ImplicitSolventKit<double> &isk,
              const NeckGeneralizedBornKit<double> &neck_gbk, double* effective_gb_radii,
              double* psi, double* sumdeijda,
              const RestraintKit<double, double2, double4_16a> &rar,
              const VirtualSiteKit<double> &vsk, const ChemicalDetailsKit &cdk,
              const ConstraintKit<double> &cnk, const StaticExclusionMaskReader &ser,
              const DynamicsControls &dyncon, int system_index = 0);
/// \}

/// \brief Carry out molecular dynamics in implicit solvent (or vacuum conditions) for a specified
///        number of steps.
///
/// \param ps                    Coordinates (positions, velocities, and forces) for all particles
///                              in the system
/// \param heat_bath             The thermostat regulating the system at a given temperature, or
///                              even regulating parts of the system at distinct temperatures
/// \param sc                    Energy tracking object
/// \param ag                    System topology (parameters for all constitutive energy terms)
/// \param neck_gbtab            "Neck" Generalized Born parameter tables
/// \param se                    Exclusion masks for non-bonded interactions among all particles
/// \param ra                    Restraint parameters for the system
/// \param dyncon                Information obtained from a &dynamics control namelist, or mocked
///                              by a developer to pass through the same input pathway
/// \param system_index          Index of the system in some larger collection of systems (this is
///                              for accessing the proper slots in the energy tracking object sc)
/// \param trajectory_file_name  Name of the trajectory file to write during production
/// \{
void dynamics(PhaseSpace *ps, Thermostat *heat_bath, ScoreCard *sc, const AtomGraph *ag,
              const NeckGeneralizedBornTable *neck_gbtab, const StaticExclusionMask *se,
              const RestraintApparatus *ra, const DynamicsControls &dyncon, int system_index = 0,
              const std::string &trajectory_file_name = std::string(""));

void dynamics(PhaseSpace *ps, Thermostat *heat_bath, ScoreCard *sc, const AtomGraph &ag,
              const NeckGeneralizedBornTable &neck_gbtab, const StaticExclusionMask &se,
              const RestraintApparatus &ra, const DynamicsControls &dyncon, int system_index = 0,
              const std::string &trajectory_file_name = std::string(""));
  
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void dynamics(PhaseSpaceSynthesis *poly_ps, CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> *cg,
              PMIGrid *pmig, ConvolutionManager *cvol, ScoreCard *sc, Thermostat *heat_bath,
              const AtomGraphSynthesis &poly_ag, const LocalExclusionMask &lem,
              const DynamicsControls &dyncon, const PrecisionControls &preccon,
              const PPPMControls &pmecon);
/// \}
  
/// \brief Carry out molecular dynamics in explicit solvent for a specified number of steps.
  
} // namespace mm
} // namespace stormm

#include "dynamics.tpp"

#endif
