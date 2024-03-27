// -*-c++-*-
#ifndef STORMM_NONBONDED_POTENTIAL_H
#define STORMM_NONBONDED_POTENTIAL_H

#include "copyright.h"
#include "Constants/generalized_born.h"
#include "DataTypes/common_types.h"
#include "Math/rounding.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"
#include "energy_abstracts.h"
#include "energy_enumerators.h"
#include "scorecard.h"
#include "soft_core_potentials.h"
#include "static_exclusionmask.h"

namespace stormm {
namespace energy {

using data_types::isSignedIntegralScalarType;
using stmath::roundUp;
using topology::AtomGraph;
using topology::ImplicitSolventKit;
using topology::NonbondedKit;
using topology::UnitCellType;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateFrameWriter;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;
using trajectory::CoordinateSeriesWriter;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;
using namespace generalized_born_defaults;

/// \brief Evaluate the non-bonded electrostatic energy using an all-to-all approach with no
///        imaging considerations (isolated boundary conditions).
///
/// Overloaded:
///   - Operate on raw pointers for the coordinates, box transformations, and forces
///   - Pass the original topology by pointer or by reference
///   - Pass the non-bonded parameter kit by value
///   - Evaluate based on a PhaseSpace object, with the option to compute forces
///   - Evaluate energy only, based on a CoordinateFrame abstract
///
/// \param ag              The system topology
/// \param nbk             Non-bonded parameters abstract taken from the original topology
/// \param se              Exclusion mask providing bits for all atom pairs
/// \param ser             Abstract of the exclusion mask providing bits for all atom pairs
/// \param ps              Coordinates and force accumulators (modified by this function)
/// \param psw             Coordinates and force accumulators (modified by this function)
/// \param cfr             Coordinates of all particles, plus box dimensions (if needed)
/// \param cfw             Coordinates of all particles, plus box dimensions (if needed)
/// \param csr             Coordinates of many frames of a particular system.  If a series is
///                        provided, the energy and state variable tracking object should be
///                        allocated to store results for all frames and the system index will be
///                        taken as the frame to evaluate.
/// \param xcrd            Cartesian X coordinates of all particles
/// \param ycrd            Cartesian Y coordinates of all particles
/// \param zcrd            Cartesian Z coordinates of all particles
/// \param xfrc            Cartesian X forces acting on all particles
/// \param yfrc            Cartesian X forces acting on all particles
/// \param zfrc            Cartesian X forces acting on all particles
/// \param umat            Box space transformation matrix
/// \param invu            Inverse transformation matrix, fractional coordinates back to real space
/// \param unit_cell       The unit cell type, i.e. triclinic
/// \param ecard           Energy components and other state variables (volume, temperature, etc.)
///                        (modified by this function)
/// \param eval_force      Flag to have forces also evaluated
/// \param system_index    Index of the system to which this energy contributes
/// \param clash_distance  Minimum distance by which any two particles must be separated in order
///                        to avoid being declared to clash, if clash detection is in effect
/// \param clash_ratio     Minimum ratio of the inter-particle distance and the pairwise van-der
///                        Waals sigma parameter by which two particles must be separated in order
///                        to avoid being declared to clash, if clash detection is in effect
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc>
double2 evaluateNonbondedEnergy(const NonbondedKit<Tcalc> nbk, const StaticExclusionMaskReader ser,
                                const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                                const double* umat, const double* invu, UnitCellType unit_cell,
                                Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, ScoreCard *ecard,
                                EvaluateForce eval_elec_force, EvaluateForce eval_vdw_force,
                                int system_index, Tcalc inv_gpos_factor = 1.0,
                                Tcalc force_factor = 1.0, Tcalc clash_distance = 0.0,
                                Tcalc clash_ratio = 0.0);

double2 evaluateNonbondedEnergy(const NonbondedKit<double> nbk,
                                const StaticExclusionMaskReader ser, PhaseSpaceWriter psw,
                                ScoreCard *ecard,
                                EvaluateForce eval_elec_force = EvaluateForce::NO,
                                EvaluateForce eval_vdw_force  = EvaluateForce::NO,
                                int system_index = 0, double clash_distance = 0.0,
				double clash_ratio = 0.0);

double2 evaluateNonbondedEnergy(const AtomGraph &ag, const StaticExclusionMask &se, PhaseSpace *ps,
                                ScoreCard *ecard,
                                EvaluateForce eval_elec_force = EvaluateForce::NO,
                                EvaluateForce eval_vdw_force  = EvaluateForce::NO,
                                int system_index = 0, double clash_distance = 0.0,
				double clash_ratio = 0.0);

double2 evaluateNonbondedEnergy(const AtomGraph *ag, const StaticExclusionMask &se, PhaseSpace *ps,
                                ScoreCard *ecard,
                                EvaluateForce eval_elec_force = EvaluateForce::NO,
                                EvaluateForce eval_vdw_force  = EvaluateForce::NO,
                                int system_index = 0, double clash_distance = 0.0,
				double clash_ratio = 0.0);

double2 evaluateNonbondedEnergy(const NonbondedKit<double> nbk,
                                const StaticExclusionMaskReader ser, CoordinateFrameReader cfr,
                                ScoreCard *ecard, int system_index = 0,
                                double clash_distance = 0.0, double clash_ratio = 0.0);

double2 evaluateNonbondedEnergy(const NonbondedKit<double> nbk,
                                const StaticExclusionMaskReader ser,
                                const CoordinateFrameWriter &cfw, ScoreCard *ecard,
                                int system_index = 0, double clash_distance = 0.0,
				double clash_ratio = 0.0);

double2 evaluateNonbondedEnergy(const AtomGraph &ag, const StaticExclusionMask &se,
                                const CoordinateFrame &cf, ScoreCard *ecard, int system_index = 0,
                                double clash_distance = 0.0, double clash_ratio = 0.0);

double2 evaluateNonbondedEnergy(const AtomGraph *ag, const StaticExclusionMask &se,
                                const CoordinateFrame &cf, ScoreCard *ecard, int system_index = 0,
                                double clash_distance = 0.0, double clash_ratio = 0.0);

template <typename Tcoord, typename Tcalc>
double2 evaluateNonbondedEnergy(const NonbondedKit<Tcalc> &nbk,
                                const StaticExclusionMaskReader &ser,
                                const CoordinateSeriesReader<Tcoord> csr, ScoreCard *ecard,
                                int system_index = 0, Tcalc clash_distance = 0.0,
				Tcalc clash_ratio = 0.0);
/// \}

/// \brief Evaluate the non-bonded Generalized Born energy of a system of particles (no cutoff is
///        applied in computation of the radii or evaluation of the energy and forces)
///
/// Overloaded:
///   - Operate on raw pointers for the coordinates, box transformations, and forces
///   - Pass the original topology by pointer or by reference,
///   - Pass the non-bonded parameter kit by value
///   - Evaluate based on a PhaseSpace object, with the option to compute forces
///   - Evaluate energy only, based on a CoordinateFrame abstract
///
/// \param ag            System topology
/// \param nbk           Non-bonded parameters abstract taken from the original topology
/// \param se            Exclusion mask providing bits for all atom pairs
/// \param ser           Abstract of the exclusion mask providing bits for all atom pairs
/// \param ngb_tables    "Neck" Generalized Born tables from some pre-loaded cache of constants
/// \param ps            Coordinates and force accumulators (modified by this function)
/// \param psw           Coordinates and force accumulators (modified by this function)
/// \param cfr           Coordinates of all particles, plus box dimensions (if needed)
/// \param cfw           Coordinates of all particles, plus box dimensions (if needed)
/// \param csr           Coordinates of many frames of a particular system.  If a series is
///                      provided, the energy and state variable tracking object should be
///                      allocated to store results for all frames and the system index will be
///                      taken as the frame to evaluate.
/// \param xcrd          Cartesian X coordinates of all particles
/// \param ycrd          Cartesian Y coordinates of all particles
/// \param zcrd          Cartesian Z coordinates of all particles
/// \param xfrc          Cartesian X forces acting on all particles
/// \param yfrc          Cartesian X forces acting on all particles
/// \param zfrc          Cartesian X forces acting on all particles
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc>
double evaluateGeneralizedBornEnergy(const NonbondedKit<Tcalc> nbk,
                                     const StaticExclusionMaskReader ser,
                                     const ImplicitSolventKit<Tcalc> isk,
                                     const NeckGeneralizedBornKit<Tcalc> ngb_kit,
                                     const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                                     Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
                                     Tforce *effective_gb_radii, Tforce *psi, Tforce *sumdeijda,
                                     ScoreCard *ecard, EvaluateForce eval_force,
                                     int system_index = 0, Tcalc inv_gpos_factor = 1.0,
                                     Tcalc force_factor = 1.0);

double evaluateGeneralizedBornEnergy(const NonbondedKit<double> nbk,
                                     const StaticExclusionMaskReader ser,
                                     const ImplicitSolventKit<double> isk,
                                     const NeckGeneralizedBornKit<double> ngb_kit,
                                     PhaseSpaceWriter psw, ScoreCard *ecard,
                                     EvaluateForce eval_force = EvaluateForce::NO,
                                     int system_index = 0);

double evaluateGeneralizedBornEnergy(const AtomGraph &ag, const StaticExclusionMask &se,
                                     const NeckGeneralizedBornTable &ngb_tables, PhaseSpace *ps,
                                     ScoreCard *ecard,
                                     EvaluateForce eval_force = EvaluateForce::NO,
                                     int system_index = 0);

double evaluateGeneralizedBornEnergy(const AtomGraph *ag, const StaticExclusionMask &se,
                                     const NeckGeneralizedBornTable &ngb_tables, PhaseSpace *ps,
                                     ScoreCard *ecard,
                                     EvaluateForce eval_force = EvaluateForce::NO,
                                     int system_index = 0);

double evaluateGeneralizedBornEnergy(const NonbondedKit<double> nbk,
                                     const StaticExclusionMaskReader ser,
                                     const ImplicitSolventKit<double> isk,
                                     const NeckGeneralizedBornKit<double> ngb_kit,
                                     const CoordinateFrameReader cfr, ScoreCard *ecard,
                                     int system_index = 0);

double evaluateGeneralizedBornEnergy(const NonbondedKit<double> nbk,
                                     const StaticExclusionMaskReader ser,
                                     const ImplicitSolventKit<double> isk,
                                     const NeckGeneralizedBornKit<double> ngb_kit,
                                     const CoordinateFrameWriter &cfw, ScoreCard *ecard,
                                     int system_index = 0);

double evaluateGeneralizedBornEnergy(const AtomGraph &ag, const StaticExclusionMask &se,
                                     const NeckGeneralizedBornTable &ngb_tables,
                                     const CoordinateFrame &cf, ScoreCard *ecard,
                                     int system_index = 0);

double evaluateGeneralizedBornEnergy(const AtomGraph *ag, const StaticExclusionMask &se,
                                     const NeckGeneralizedBornTable &ngb_tables,
                                     const CoordinateFrame &cf, ScoreCard *ecard,
                                     int system_index = 0);

template <typename Tcoord, typename Tcalc>
double evaluateGeneralizedBornEnergy(const NonbondedKit<Tcalc> nbk,
                                     const StaticExclusionMaskReader ser,
                                     const ImplicitSolventKit<Tcalc> isk,
                                     const NeckGeneralizedBornKit<Tcalc> ngb_kit,
                                     const CoordinateSeriesReader<Tcoord> csr, ScoreCard *ecard,
                                     int system_index = 0, int force_scale_bits = 23);
/// \}

} // namespace energy
} // namespace stormm

#include "nonbonded_potential.tpp"

#endif
