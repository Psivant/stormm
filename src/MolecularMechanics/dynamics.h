// -*-c++-*-
#ifndef STORMM_DYNAMICS_LIBRARY_H
#define STORMM_DYNAMICS_LIBRARY_H

#include "copyright.h"
#include "Constants/generalized_born.h"
#include "DataTypes/stormm_vector_types.h"
#include "MolecularMechanics/mm_evaluation.h"
#include "Potential/static_exclusionmask.h"
#include "Restraints/restraint_apparatus.h"
#include "Structure/rattle.h"
#include "Structure/virtual_site_handling.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/integration.h"
#include "Trajectory/phasespace.h"
#include "Trajectory/thermostat.h"
#include "kinetic.h"

namespace stormm {
namespace mm {

using energy::ScoreCard;
using energy::StaticExclusionMask;
using energy::StaticExclusionMaskReader;
using namelist::DynamicsControls;
using restraints::RestraintApparatus;
using restraints::RestraintKit;
using structure::placeVirtualSites;
using structure::rattlePositions;
using structure::rattleVelocities;
using structure::transmitVirtualSiteForces;
using topology::AtomGraph;
using topology::ConstraintKit;
using topology::ImplicitSolventKit;
using topology::NonbondedKit;
using topology::ValenceKit;
using topology::VirtualSiteKit;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;
using trajectory::Thermostat;
using trajectory::ThermostatWriter;
using trajectory::velocityVerletVelocityUpdate;
using trajectory::velocityVerletCoordinateUpdate;
using namespace generalized_born_defaults;
  
/// \brief Run the canonical MD dynamics step for systems in implicit solvent.
///
/// \param xcrd
/// \param ycrd
/// \param zcrd
/// \param xvel
/// \param yvel
/// \param zvel
/// \param xfrc
/// \param yfrc
/// \param zfrc
/// \param xalt
/// \param yalt
/// \param zalt
/// \param vxalt
/// \param vyalt
/// \param vzalt
/// \param fxalt
/// \param fyalt
/// \param fzalt
/// \param psw
/// \param tstr
/// \param tst
/// \param vk
/// \param nbk
/// \param isk
/// \param neck_gbk
/// \param effective_gb_radii  Array to accumulate the effective GB radii for all atoms
/// \param psi                 Array to hold the GB Psi values for all atoms
/// \param sumdeijda           Array to accumulate the GB derivative factors for all atoms
/// \param rar   
/// \param vsk
/// \param cnk
/// \param ser
/// \param dyncon
/// \param nrg_scale_bits
/// \param gpos_scale_factor   Scaling factor on 
/// \param vel_scale_factor
/// \param frc_scale_factor    
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

void dynaStep(PhaseSpaceWriter *psw, ScoreCard *sc, const ThermostatWriter<double> &tstr,
              const ValenceKit<double> &vk, const NonbondedKit<double> &nbk,
              const ImplicitSolventKit<double> &isk,
              const NeckGeneralizedBornKit<double> &neck_gbk, double* effective_gb_radii,
              double* psi, double* sumdeijda, const RestraintKit<double, double2, double4> &rar,
              const VirtualSiteKit<double> &vsk, const ChemicalDetailsKit &cdk,
              const ConstraintKit<double> &cnk, const StaticExclusionMaskReader &ser,
              const DynamicsControls &dyncon, int system_index = 0);
/// \}

/// \brief Carry out molecular dynamics in implicit solvent (or vacuum conditions) for a specified
///        number of steps.
///
/// \param ps
/// \param heat_bath
/// \param sc
/// \param ag
/// \param neck_gbtab
/// \param se
/// \param ra
/// \param dyncon
/// \param system_index
/// \{
void dynamics(PhaseSpace *ps, Thermostat *heat_bath, ScoreCard *sc, const AtomGraph *ag,
              const NeckGeneralizedBornTable *neck_gbtab, const StaticExclusionMask *se,
              const RestraintApparatus *ra, const DynamicsControls &dyncon, int system_index = 0,
              const std::string &trajectory_file_name = std::string(""),
              const std::string &restart_file_name = std::string(""));

void dynamics(PhaseSpace *ps, Thermostat *heat_bath, ScoreCard *sc, const AtomGraph &ag,
              const NeckGeneralizedBornTable &neck_gbtab, const StaticExclusionMask &se,
              const RestraintApparatus &ra, const DynamicsControls &dyncon, int system_index = 0,
              const std::string &trajectory_file_name = std::string(""),
              const std::string &restart_file_name = std::string(""));
/// \}
  
} // namespace mm
} // namespace stormm

#include "dynamics.tpp"

#endif
