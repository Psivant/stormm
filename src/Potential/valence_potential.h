// -*-c++-*-
#ifndef STORMM_VALENCE_POTENTIAL_H
#define STORMM_VALENCE_POTENTIAL_H

#include <cmath>
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "Math/matrix_ops.h"
#include "Math/rounding.h"
#include "Math/vector_ops.h"
#include "Restraints/restraint_apparatus.h"
#include "Restraints/restraint_util.h"
#include "Structure/local_arrangement.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"
#include "energy_enumerators.h"
#include "scorecard.h"
#include "soft_core_potentials.h"

namespace stormm {
namespace energy {

using data_types::isSignedIntegralScalarType;
using stmath::angleVerification;
using stmath::crossProduct;
using stmath::matrixMultiply;
using stmath::matrixVectorMultiply;
using stmath::roundUp;
using restraints::computeRestraintMixture;
using restraints::RestraintApparatus;
using restraints::RestraintKit;
using restraints::restraintDelta;
using structure::imageCoordinates;
using structure::ImagingMethod;
using symbols::asymptotic_to_one_f;
using symbols::asymptotic_to_one_lf;
using symbols::inverse_twopi;
using symbols::inverse_twopi_f;
using symbols::inverse_one_minus_asymptote_f;
using symbols::inverse_one_minus_asymptote_lf;
using symbols::pi;
using symbols::pi_f;
using symbols::twopi;
using symbols::twopi_f;
using topology::AtomGraph;
using topology::ValenceKit;
using topology::NonbondedKit;
using topology::TorsionKind;
using topology::UnitCellType;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateFrameWriter;
using trajectory::CoordinateSeriesReader;
using trajectory::CoordinateSeriesWriter;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;

/// \brief Evaluate the energy and forces dur to a harmonic stretching term (both harmonic bonds
///        and Urey-Bradley terms can be evaluated with this function).
///
/// \param i_atom           Index of the first atom in the interaction
/// \param j_atom           Index of the second atom in the interaction
/// \param stiffness        Stiffness of the stretching potential (units of kcal/mol-Angstrom)
/// \param equilibrium      Equilibrium value of the stretching function (units of Angstroms)
/// \param xcrd             Cartesian X coordinates of all particles
/// \param ycrd             Cartesian Y coordinates of all particles
/// \param zcrd             Cartesian Z coordinates of all particles
/// \param umat             Box space transformation matrix, real coordinates to fractional space
/// \param invu             Inverse transformation matrix, fractional coordinates to real space
/// \param unit_cell        The unit cell type, i.e. triclinic
/// \param xfrc             Cartesian X forces acting on all particles
/// \param yfrc             Cartesian Y forces acting on all particles
/// \param zfrc             Cartesian Z forces acting on all particles
/// \param eval_force       Flag to have forces also evaluated
/// \param inv_gpos_factor  Inverse positional scaling factor (for signed integer, fixed-precision
///                         coordinate representations)
/// \param force_factor     Force scaling factor for fixed-precision force accumulation
template <typename Tcoord, typename Tforce, typename Tcalc>
Tcalc evalHarmonicStretch(int i_atom, int j_atom, Tcalc stiffness, Tcalc equilibrium,
                          const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                          const double* umat, const double* invu, UnitCellType unit_cell,
                          Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, EvaluateForce eval_force,
                          Tcalc inv_gpos_factor = 1.0, Tcalc force_factor = 1.0);

/// \brief Evaluate the bond energy contributions based on a topology and a coordinate set.
///        These simple routines can serve as a check on much more complex routines involving
///        streamlined data structures and GPU execution.
///
/// Overloaded:
///   - Evaluate based on raw pointers to coordinates, box transformations, and forces
///   - Evaluate based on a PhaseSpace object, with the option to compute and store forces
///   - Evaluate energy only based on a CoordinateFrame or CoordinateSeries abstract
///   - Pass a topology by pointer, by reference, or just the ValenceKit abstract by value
///
/// \param ag            System topology
/// \param vk            Valence parameters abstract from the system topology
/// \param ps            Coordinates, box size, and force accumulators (modified by this function)
/// \param psw           Coordinates, box size, and force accumulators (modified by this function)
/// \param cfr           Coordinates of all particles, plus box dimensions (if needed)
/// \param cfw           Coordinates of all particles, plus box dimensions (if needed)
/// \param csr           Coordinates of a series of snapshots, each containing all particles.  If
///                      using such a series the energy tracking object should be allocated to
///                      store results for all frames and the system index will be taken as the
///                      frame number of interest.
/// \param xcrd          Cartesian X coordinates of all particles
/// \param ycrd          Cartesian Y coordinates of all particles
/// \param zcrd          Cartesian Z coordinates of all particles
/// \param xfrc          Cartesian X forces acting on all particles
/// \param yfrc          Cartesian X forces acting on all particles
/// \param zfrc          Cartesian X forces acting on all particles
/// \param umat          Box space transformation matrix
/// \param invu          Inverse transformation matrix, fractional coordinates back to real space
/// \param unit_cell     The unit cell type, i.e. triclinic
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc>
double evaluateBondTerms(const ValenceKit<Tcalc> vk, const Tcoord* xcrd, const Tcoord* ycrd,
                         const Tcoord* zcrd, const double* umat, const double* invu,
                         UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
                         ScoreCard *ecard, EvaluateForce eval_force = EvaluateForce::NO,
                         int system_index = 0, Tcalc inv_gpos_factor = 1.0,
                         Tcalc force_factor = 1.0);
  
double evaluateBondTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                         EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateBondTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                         EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateBondTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                         EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateBondTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                         ScoreCard *ecard, int system_index = 0);

double evaluateBondTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                         ScoreCard *ecard, int system_index = 0);

double evaluateBondTerms(const AtomGraph &ag, const CoordinateFrame &cf,
                         ScoreCard *ecard, int system_index = 0);

double evaluateBondTerms(const AtomGraph *ag, const CoordinateFrame &cf,
                         ScoreCard *ecard, int system_index = 0);

template <typename Tcoord, typename Tcalc>
double evaluateBondTerms(const ValenceKit<Tcalc> vk, const CoordinateSeriesReader<Tcoord> csr,
                         ScoreCard *ecard, int system_index = 0, int force_scale_bits = 23);

template <typename Tcoord, typename Tcalc>
double evaluateBondTerms(const ValenceKit<Tcalc> vk, const CoordinateSeriesWriter<Tcoord> csr,
                         ScoreCard *ecard, int system_index = 0, int force_scale_bits = 23);
/// \}

/// \brief Evaluate the energy and forces due to a harmonic bending interaction.  Parameters for
///        this function follow evalHarmonicStretch, with the addition of:
///
/// \param k_atom       Index of the third atom in the interaction
/// \param stiffness    Stiffness constant for the bending function, units of kcal/mol-radian
/// \param equilibrium  Equilibrium value for the bending function, units of radians
template <typename Tcoord, typename Tforce, typename Tcalc>
Tcalc evalHarmonicBend(int i_atom, int j_atom, int k_atom, Tcalc stiffness, Tcalc equilibrium,
                       const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                       const double* umat, const double* invu, UnitCellType unit_cell,
                       Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, EvaluateForce eval_force,
                       Tcalc inv_gpos_factor = 1.0, Tcalc force_factor = 1.0);

/// \brief Evaluate the angle bending energy contributions with a simple routine based on a
///        topology and a coordinate set.  These simple routines can serve as a check on much more
///        complex routines involving streamlined data structures and GPU execution.
///
/// Overloaded:
///   - Evaluate based on raw pointers to coordinates, box transformations, and forces
///   - Evaluate based on a PhaseSpace object, with the option to compute and store forces
///   - Evaluate energy only based on a CoordinateFrame or CoordinateSeries abstract
///   - Pass a topology by pointer, by reference, or just the ValenceKit abstract by value
///
/// \param ag            System topology
/// \param vk            Valence parameters abstract from the system topology
/// \param ps            Coordinates, box size, and force accumulators (modified by this function)
/// \param psw           Coordinates, box size, and force accumulators (modified by this function)
/// \param cfr           Coordinates of all particles, plus box dimensions (if needed)
/// \param cfw           Coordinates of all particles, plus box dimensions (if needed)
/// \param csr           Coordinates of a series of snapshots, each containing all particles.  If
///                      using such a series the energy tracking object should be allocated to
///                      store results for all frames and the system index will be taken as the
///                      frame number of interest.
/// \param xcrd          Cartesian X coordinates of all particles
/// \param ycrd          Cartesian Y coordinates of all particles
/// \param zcrd          Cartesian Z coordinates of all particles
/// \param xfrc          Cartesian X forces acting on all particles
/// \param yfrc          Cartesian X forces acting on all particles
/// \param zfrc          Cartesian X forces acting on all particles
/// \param umat          Box space transformation matrix
/// \param invu          Inverse transformation matrix, fractional coordinates back to real space
/// \param unit_cell     The unit cell type, i.e. triclinic
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc>
double evaluateAngleTerms(const ValenceKit<Tcalc> vk, const Tcoord* xcrd, const Tcoord* ycrd,
                          const Tcoord* zcrd, const double* umat, const double* invu,
                          UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
                          ScoreCard *ecard, EvaluateForce eval_force = EvaluateForce::NO,
                          int system_index = 0, Tcalc inv_gpos_factor = 1.0,
                          Tcalc force_factor = 1.0);

double evaluateAngleTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                          EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);
                          
double evaluateAngleTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                          EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateAngleTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                          EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateAngleTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                          ScoreCard *ecard, int system_index = 0);

double evaluateAngleTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                          ScoreCard *ecard, int system_index = 0);

double evaluateAngleTerms(const AtomGraph &ag, const CoordinateFrame &cf, ScoreCard *ecard,
                          int system_index = 0);

double evaluateAngleTerms(const AtomGraph *ag, const CoordinateFrame &cf, ScoreCard *ecard,
                          int system_index = 0);

template <typename Tcoord, typename Tcalc>
double evaluateAngleTerms(const ValenceKit<Tcalc> vk, const CoordinateSeriesReader<Tcoord> csr,
                          ScoreCard *ecard, int system_index = 0, int force_scale_bits = 23);

template <typename Tcoord, typename Tcalc>
double evaluateAngleTerms(const ValenceKit<Tcalc> vk, const CoordinateSeriesWriter<Tcoord> csr,
                          ScoreCard *ecard, int system_index = 0, int force_scale_bits = 23);
/// \}
  
/// \brief Evalaute the energy and forces due to a cosine-based or harmonic dihedral term.
///        Parameters for this function follow evalHarmonicBend, with the addition of:
///
/// \param l_atom       Index of the fourth atom in the interaction
/// \param amplitude    Amplitude of the cosine wave, units of kcal/mol
/// \param phase_angle  Phase shift of the cosine function argument, units of radians
/// \param frequency    Periodicity (frequency) in the cosine function argument
/// \param kind         Specifies that the dihedral potential is either a cosine series term or a
///                     harmonic penalty function
template <typename Tcoord, typename Tforce, typename Tcalc>
Tcalc evalDihedralTwist(int i_atom, int j_atom, int k_atom, int l_atom, Tcalc amplitude,
                        Tcalc phase_angle, Tcalc frequency, DihedralStyle kind,
                        const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                        const double* umat, const double* invu, UnitCellType unit_cell,
                        Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, EvaluateForce eval_force,
                        Tcalc inv_gpos_factor = 1.0, Tcalc force_factor = 1.0);

/// \brief Evaluate the proper and improper dihedral energy contributions with a simple routine
///        based on a topology and a PhaseSpace object to store forces in double precision.  The
///        contributions of proper and improper dihedrals are stored in the x and y components
///        of the resulting double-precision tuple, respectively.  These results can be compared
///        to the fixed-precision accumulators in the energy tracking object.
///
/// Overloaded:
///   - Evaluate based on raw pointers to coordinates, box transformations, and forces
///   - Evaluate based on a PhaseSpace object, with the option to compute and store forces
///   - Evaluate energy only based on a CoordinateFrame abstract
///   - Pass a topology by pointer, by reference, or just the ValenceKit abstract by value
///
/// \param ag            System topology
/// \param vk            Valence parameters abstract from the system topology
/// \param ps            Coordinates, box size, and force accumulators (modified by this function)
/// \param psw           Coordinates, box size, and force accumulators (modified by this function)
/// \param cfr           Coordinates of all particles, plus box dimensions (if needed)
/// \param cfw           Coordinates of all particles, plus box dimensions (if needed)
/// \param xcrd          Cartesian X coordinates of all particles
/// \param ycrd          Cartesian Y coordinates of all particles
/// \param zcrd          Cartesian Z coordinates of all particles
/// \param xfrc          Cartesian X forces acting on all particles
/// \param yfrc          Cartesian X forces acting on all particles
/// \param zfrc          Cartesian X forces acting on all particles
/// \param umat          Box space transformation matrix
/// \param invu          Inverse transformation matrix, fractional coordinates back to real space
/// \param unit_cell     The unit cell type, i.e. triclinic
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc>
double2 evaluateDihedralTerms(const ValenceKit<Tcalc> vk, const Tcoord* xcrd, const Tcoord* ycrd,
                              const Tcoord* zcrd, const double* umat, const double* invu,
                              UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
                              ScoreCard *ecard, EvaluateForce eval_force = EvaluateForce::NO,
                              int system_index = 0, Tcalc inv_gpos_factor = 1.0,
                              Tcalc force_factor = 1.0);

double2 evaluateDihedralTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                              EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double2 evaluateDihedralTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                              EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double2 evaluateDihedralTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                              EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double2 evaluateDihedralTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                              ScoreCard *ecard, int system_index = 0);

double2 evaluateDihedralTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                              ScoreCard *ecard, int system_index = 0);

double2 evaluateDihedralTerms(const AtomGraph &ag, const CoordinateFrame &cf, ScoreCard *ecard,
                              int system_index = 0);

double2 evaluateDihedralTerms(const AtomGraph *ag, const CoordinateFrame &cf, ScoreCard *ecard,
                              int system_index = 0);

template <typename Tcoord, typename Tcalc>
double2 evaluateDihedralTerms(const ValenceKit<Tcalc> vk, const CoordinateSeriesReader<Tcoord> csr,
                              ScoreCard *ecard, int system_index = 0, int force_scale_bits = 23);

template <typename Tcoord, typename Tcalc>
double2 evaluateDihedralTerms(const ValenceKit<Tcalc> vk, const CoordinateSeriesWriter<Tcoord> csr,
                              ScoreCard *ecard, int system_index = 0, int force_scale_bits = 23);
/// \}
  
/// \brief Evaluate Urey-Bradley harmonic angle interactions with a simple routine.  This looks
///        almost exactly like the bond computations but it kept in this separate routine for
///        simplicity of bookkeeping.  The resulting energy accumulated in double-precision can be
///        compared to the result obtained with fixed-precision accumulation.
///
/// Overloaded:
///   - Evaluate based on raw pointers to coordinates, box transformations, and forces
///   - Evaluate based on a PhaseSpace object, with the option to compute and store forces
///   - Evaluate energy only based on a CoordinateFrame abstract
///   - Pass a topology by pointer, by reference, or just the ValenceKit abstract by value
///
/// \param ag            System topology
/// \param vk            Valence parameters abstract from the system topology
/// \param ps            Coordinates, box size, and force accumulators (modified by this function)
/// \param psw           Coordinates, box size, and force accumulators (modified by this function)
/// \param cfr           Coordinates of all particles, plus box dimensions (if needed)
/// \param cfw           Coordinates of all particles, plus box dimensions (if needed)
/// \param xcrd          Cartesian X coordinates of all particles
/// \param ycrd          Cartesian Y coordinates of all particles
/// \param zcrd          Cartesian Z coordinates of all particles
/// \param xfrc          Cartesian X forces acting on all particles
/// \param yfrc          Cartesian X forces acting on all particles
/// \param zfrc          Cartesian X forces acting on all particles
/// \param umat          Box space transformation matrix
/// \param invu          Inverse transformation matrix, fractional coordinates back to real space
/// \param unit_cell     The unit cell type, i.e. triclinic
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc>
double evaluateUreyBradleyTerms(const ValenceKit<Tcalc> vk, const Tcoord* xcrd,
                                const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                                const double* invu, const UnitCellType unit_cell, Tforce* xfrc,
                                Tforce* yfrc, Tforce* zfrc, ScoreCard *ecard,
                                EvaluateForce eval_force = EvaluateForce::NO,
                                int system_index = 0, Tcalc inv_gpos_factor = 1.0,
                                Tcalc force_factor = 1.0);

double evaluateUreyBradleyTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw,
                                ScoreCard *ecard, EvaluateForce eval_force = EvaluateForce::NO,
                                int system_index = 0);

double evaluateUreyBradleyTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                                EvaluateForce eval_force = EvaluateForce::NO,
                                int system_index = 0);

double evaluateUreyBradleyTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                                EvaluateForce eval_force = EvaluateForce::NO,
                                int system_index = 0);

double evaluateUreyBradleyTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                                ScoreCard *ecard, int system_index = 0);

double evaluateUreyBradleyTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                                ScoreCard *ecard, int system_index = 0);

double evaluateUreyBradleyTerms(const AtomGraph &ag, const CoordinateFrame &cf, ScoreCard *ecard,
                                int system_index = 0);

double evaluateUreyBradleyTerms(const AtomGraph *ag, const CoordinateFrame &cf, ScoreCard *ecard,
                                int system_index = 0);

template <typename Tcoord, typename Tcalc>
double evaluateUreyBradleyTerms(const ValenceKit<Tcalc> vk,
                                const CoordinateSeriesReader<Tcoord> csr, ScoreCard *ecard,
                                int system_index = 0, int force_scale_bits = 23);

template <typename Tcoord, typename Tcalc>
double evaluateUreyBradleyTerms(const ValenceKit<Tcalc> vk,
                                const CoordinateSeriesWriter<Tcoord> csr, ScoreCard *ecard,
                                int system_index = 0, int force_scale_bits = 23);
/// \}
  
/// \brief Evaluate CHARMM harmonic improper dihedral terms with a simple routine.  This
///        contributes to a separate improper energy accumulator, CHARMM_IMPROPER, in the
///        ScoreCard object.  The double-precision output of this function can be compared to the
///        result in fixed-precision accumulation.
///
/// Overloaded:
///   - Evaluate based on raw pointers to coordinates, box transformations, and forces
///   - Evaluate based on a PhaseSpace object, with the option to compute and store forces
///   - Evaluate energy only based on a CoordinateFrame abstract
///   - Pass a topology by pointer, by reference, or just the ValenceKit abstract by value
///
/// \param ag            System topology
/// \param vk            Valence parameters abstract from the system topology
/// \param ps            Coordinates, box size, and force accumulators (modified by this function)
/// \param psw           Coordinates, box size, and force accumulators (modified by this function)
/// \param cfr           Coordinates of all particles, plus box dimensions (if needed)
/// \param cfw           Coordinates of all particles, plus box dimensions (if needed)
/// \param xcrd          Cartesian X coordinates of all particles
/// \param ycrd          Cartesian Y coordinates of all particles
/// \param zcrd          Cartesian Z coordinates of all particles
/// \param xfrc          Cartesian X forces acting on all particles
/// \param yfrc          Cartesian X forces acting on all particles
/// \param zfrc          Cartesian X forces acting on all particles
/// \param umat          Box space transformation matrix
/// \param invu          Inverse transformation matrix, fractional coordinates back to real space
/// \param unit_cell     The unit cell type, i.e. triclinic
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc>
double evaluateCharmmImproperTerms(const ValenceKit<Tcalc> vk, const Tcoord* xcrd,
                                   const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                                   const double* invu, UnitCellType unit_cell, Tforce* xfrc,
                                   Tforce* yfrc, Tforce* zfrc, ScoreCard *ecard,
                                   EvaluateForce eval_force = EvaluateForce::NO,
                                   int system_index = 0, Tcalc inv_gpos_factor = 1.0,
                                   Tcalc force_factor = 1.0);

double evaluateCharmmImproperTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw,
                                   ScoreCard *ecard, EvaluateForce eval_force = EvaluateForce::NO,
                                   int system_index = 0);

double evaluateCharmmImproperTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                                   EvaluateForce eval_force = EvaluateForce::NO,
                                   int system_index = 0);

double evaluateCharmmImproperTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                                   EvaluateForce eval_force = EvaluateForce::NO,
                                   int system_index = 0);

double evaluateCharmmImproperTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                                   ScoreCard *ecard, int system_index = 0);

double evaluateCharmmImproperTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                                   ScoreCard *ecard, int system_index = 0);

double evaluateCharmmImproperTerms(const AtomGraph &ag, const CoordinateFrame &cf,
                                   ScoreCard *ecard, int system_index = 0);

double evaluateCharmmImproperTerms(const AtomGraph *ag, const CoordinateFrame &cf,
                                   ScoreCard *ecard, int system_index = 0);

template <typename Tcoord, typename Tcalc>
double evaluateCharmmImproperTerms(const ValenceKit<Tcalc> vk,
                                   const CoordinateSeriesReader<Tcoord> csr, ScoreCard *ecard,
                                   int system_index = 0, int force_scale_bits = 23);

template <typename Tcoord, typename Tcalc>
double evaluateCharmmImproperTerms(const ValenceKit<Tcalc> vk,
                                   const CoordinateSeriesWriter<Tcoord> csr, ScoreCard *ecard,
                                   int system_index = 0, int force_scale_bits = 23);
/// \}

/// \brief Evaluate a single CHARMM CMAP two-dimensional potential and force term.  This procedure
///        is abstracted to avoid code replication.  Other lengthy procedures like the dihedral
///        computation are computed in slightly different ways that would make such encapsulation
///        difficult.
///
/// \param cmap_patches       Data array of all CMAP surfaces, in patch format (pre-computed Axy
///                           coefficients for every grid element)
/// \param cmap_patch_bounds  Bounds array for cmap_patches, with demarcations for each individual
///                           CMAP (not each individual patch)
/// \param surf_idx           The relevant surface, an index into cmap_patch_bounds
/// \param surf_dim           Dimension of the CMAP surface of interest (periodicity with 2 x pi
///                           extent along each axis is assumed)
/// \param i_atom             Atom I in the interaction (this still assumes that there are two
///                           dihedral terms with three overlapping atoms)
/// \param j_atom             Atom J in the interaction
/// \param k_atom             Atom K in the interaction
/// \param l_atom             Atom L in the interaction
/// \param m_atom             Atom M in the interaction
/// \param xcrd               Cartesian X coordinates of all particles
/// \param ycrd               Cartesian Y coordinates of all particles
/// \param zcrd               Cartesian Z coordinates of all particles
/// \param umat               Box space transformation matrix
/// \param invu               Inverse transformation matrix, fractional coordinates to real space
/// \param unit_cell          The unit cell type, i.e. triclinic
/// \param xfrc               Cartesian X forces acting on all particles
/// \param yfrc               Cartesian X forces acting on all particles
/// \param zfrc               Cartesian X forces acting on all particles
/// \param eval_force         Flag to have forces also evaluated
template <typename Tcoord, typename Tforce, typename Tcalc>
Tcalc evalCmap(const Tcalc* cmap_patches, const int* cmap_patch_bounds, int surf_idx,
               int surf_dim, int i_atom, int j_atom, int k_atom, int l_atom, int m_atom,
               const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
               const double* invu, UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
               Tforce* zfrc, EvaluateForce eval_force, Tcalc inv_gpos_factor = 1.0,
               Tcalc force_factor = 1.0);
  
/// \brief Evaluate CHARMM CMAP two-dimensional cubic spline potentials.  As with all other
///        functions in this library, the results of double-precision accumulation are returned
///        for comparison to fixed-precision results.
///
/// Overloaded:
///   - Evaluate based on raw pointers to coordinates, box transformations, and forces
///   - Evaluate based on a PhaseSpace object, with the option to compute and store forces
///   - Evaluate energy only based on a CoordinateFrame abstract
///   - Pass a topology by pointer, by reference, or just the ValenceKit abstract by value
///
/// \param ag            System topology
/// \param vk            Valence parameters abstract from the system topology
/// \param ps            Coordinates, box size, and force accumulators (modified by this function)
/// \param psw           Coordinates, box size, and force accumulators (modified by this function)
/// \param cfr           Coordinates of all particles, plus box dimensions (if needed)
/// \param cfw           Coordinates of all particles, plus box dimensions (if needed)
/// \param xcrd          Cartesian X coordinates of all particles
/// \param ycrd          Cartesian Y coordinates of all particles
/// \param zcrd          Cartesian Z coordinates of all particles
/// \param umat          Box space transformation matrix
/// \param invu          Inverse transformation matrix, fractional coordinates back to real space
/// \param unit_cell     The unit cell type, i.e. triclinic
/// \param xfrc          Cartesian X forces acting on all particles
/// \param yfrc          Cartesian Y forces acting on all particles
/// \param zfrc          Cartesian Z forces acting on all particles
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc>
double evaluateCmapTerms(const ValenceKit<Tcalc> vk, const Tcoord* xcrd, const Tcoord* ycrd,
                         const Tcoord* zcrd, const double* umat, const double* invu,
                         UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
                         ScoreCard *ecard, EvaluateForce eval_force = EvaluateForce::NO,
                         int system_index = 0, Tcalc inv_gpos_factor = 1.0,
                         Tcalc force_factor = 1.0);

double evaluateCmapTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                         EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateCmapTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                         EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateCmapTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                         EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0);

double evaluateCmapTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                         ScoreCard *ecard, int system_index = 0);

double evaluateCmapTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                         ScoreCard *ecard, int system_index = 0);

double evaluateCmapTerms(const AtomGraph &ag, const CoordinateFrame &cf, ScoreCard *ecard,
                         int system_index = 0);

double evaluateCmapTerms(const AtomGraph *ag, const CoordinateFrame &cf, ScoreCard *ecard,
                         int system_index = 0);

template <typename Tcoord, typename Tcalc>
double evaluateCmapTerms(const ValenceKit<Tcalc> vk, const CoordinateSeriesReader<Tcoord> csr,
                         ScoreCard *ecard, int system_index = 0, int force_scale_bits = 23);

template <typename Tcoord, typename Tcalc>
double evaluateCmapTerms(const ValenceKit<Tcalc> vk, const CoordinateSeriesWriter<Tcoord> csr,
                         ScoreCard *ecard, int system_index = 0, int force_scale_bits = 23);
/// \}

/// \brief Evaluate a 1:4 pair interaction according to a pair of attenuation factors for
///        electrostatic and van-der Waals interactions.  Both dihedral-bound and inferred
///        1:4 pairs can be evaluated with this one routine.
///
/// Overloaded:
///   - Use a Lennard-Jones parameter table offset, or do not
///
/// \param i_atom               The first atom in the pair
/// \param l_atom               The second atom in the pair
/// \param attn_idx             The attenuation index into arrays of electrostatic and van-der
///                             Waals scaling factors
/// \param coulomb_constant     Coulomb's constant for scaling electrostatic interactions
/// \param lj_param_idx         Lennard-Jones parameter index numbers
/// \param attn14_elec_factors  Unique electrostatic 1:4 scaling factors
/// \param attn14_vdw_factors   Unique van-der Waals 1:4 scaling factors
/// \param lja_14_coeff         Lennard-Jones A coefficients, in U = (A / r^12) - (B / r^6)
/// \param ljb_14_coeff         Lennard-Jones B coefficients, in U = (A / r^12) - (B / r^6)
/// \param lj_14_sigma          Lennard-Jones sigma coefficients, pre-calculated for convenience
/// \param ljtab_offset         Offset for the system's Lennard-Jones coefficient tables
/// \param n_lj_types           Number of unique Lennard-Jones types, the rank of the A and B
///                             coefficient matrices above
/// \param xcrd                 Cartesian X coordinates of all particles
/// \param ycrd                 Cartesian Y coordinates of all particles
/// \param zcrd                 Cartesian Z coordinates of all particles
/// \param xfrc                 Cartesian X forces acting on all particles
/// \param yfrc                 Cartesian Y forces acting on all particles
/// \param zfrc                 Cartesian Z forces acting on all particles
/// \param umat                 Box space transformation matrix
/// \param invu                 Inverse transformation, fractional coordinates back to real space
/// \param unit_cell            The unit cell type, i.e. triclinic
/// \param eval_elec_force      Flag to have electrostatic forces evaluated
/// \param eval_vdw_force       Flag to have van-der Waals (Lennard-Jones) forces evaluated
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc>
Vec2<Tcalc> evaluateAttenuated14Pair(int i_atom, int l_atom, int attn_idx, Tcalc coulomb_constant,
                                     const Tcalc* charges, const int* lj_param_idx,
                                     const Tcalc* attn14_elec_factors,
                                     const Tcalc* attn14_vdw_factors, const Tcalc* lja_14_coeff,
                                     const Tcalc* ljb_14_coeff, const Tcalc* lj_14_sigma,
                                     int ljtab_offset, int n_lj_types, const Tcoord* xcrd,
                                     const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                                     const double* invu, UnitCellType unit_cell, Tforce* xfrc,
                                     Tforce* yfrc, Tforce* zfrc, EvaluateForce eval_elec_force,
                                     EvaluateForce eval_vdw_force, Tcalc inv_gpos_factor = 1.0,
                                     Tcalc force_factor = 1.0, Tcalc clash_minimum_distance = 0.0,
                                     Tcalc clash_ratio = 0.0);

template <typename Tcoord, typename Tforce, typename Tcalc>
Vec2<Tcalc> evaluateAttenuated14Pair(int i_atom, int l_atom, int attn_idx, Tcalc coulomb_constant,
                                     const Tcalc* charges, const int* lj_param_idx,
                                     const Tcalc* attn14_elec_factors,
                                     const Tcalc* attn14_vdw_factors, const Tcalc* lja_14_coeff,
                                     const Tcalc* ljb_14_coeff, const Tcalc* lj_14_sigma,
                                     int n_lj_types, const Tcoord* xcrd, const Tcoord* ycrd,
                                     const Tcoord* zcrd, const double* umat, const double* invu,
                                     UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                                     Tforce* zfrc, EvaluateForce eval_elec_force,
                                     EvaluateForce eval_vdw_force, Tcalc inv_gpos_factor = 1.0,
                                     Tcalc force_factor = 1.0, Tcalc clash_minimum_distance = 0.0,
                                     Tcalc clash_ratio = 0.0);
/// \}
  
/// \brief Evaluate 1:4 non-bonded pair interactions.  This requires a suprising amount of
///        bookkeeping to make it performant, but the result is straightforward and this reference
///        routine will work from that setup.
///
/// Overloaded:
///   - Evaluate based on raw pointers to coordinates, box transformations, and forces
///   - Evaluate based on a PhaseSpace object, with the option to compute and store forces
///   - Evaluate energy only based on a CoordinateFrame abstract
///   - Pass a topology by pointer, by reference, or just the ValenceKit abstract by value
///
/// \param ag               System topology
/// \param vk               Valence parameters abstract from the system topology
/// \param nbk              Non-bonded parameters abstract from the system topology
/// \param ps               Coordinates, box size, and force accumulators (modified by this
///                         function)
/// \param psw              Coordinates, box size, and force accumulators (modified by this
///                         function)
/// \param cfr              Coordinates of all particles, plus box dimensions (if needed)
/// \param cfw              Coordinates of all particles, plus box dimensions (if needed)
/// \param xcrd             Cartesian X coordinates of all particles
/// \param ycrd             Cartesian Y coordinates of all particles
/// \param zcrd             Cartesian Z coordinates of all particles
/// \param umat             Box space transformation matrix
/// \param invu             Inverse transformation, fractional coordinates back to real space
/// \param unit_cell        The unit cell type, i.e. triclinic
/// \param xfrc             Cartesian X forces acting on all particles
/// \param yfrc             Cartesian Y forces acting on all particles
/// \param zfrc             Cartesian Z forces acting on all particles
/// \param ecard            Energy components and other state variables (volume, temperature, etc.)
///                         (modified by this function)
/// \param eval_elec_force  Flag to have electrostatic forces evaluated
/// \param eval_vdw_force   Flag to have van-der Waals (Lennard-Jones) forces evaluated
/// \param system_index     Index of the system to which this energy contributes
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc>
double2 evaluateAttenuated14Terms(const ValenceKit<Tcalc> vk, const NonbondedKit<Tcalc> nbk,
                                  const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                                  const double* umat, const double* invu,
                                  UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                                  Tforce* zfrc, ScoreCard *ecard,
                                  EvaluateForce eval_elec_force = EvaluateForce::NO,
                                  EvaluateForce eval_vdw_force = EvaluateForce::NO,
                                  int system_index = 0, Tcalc inv_gpos_factor = 1.0,
                                  Tcalc force_factor = 1.0, Tcalc clash_minimum_distance = 0.0,
                                  Tcalc clash_ratio = 0.0);

double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  PhaseSpaceWriter psw, ScoreCard *ecard,
                                  EvaluateForce eval_elec_force = EvaluateForce::NO,
                                  EvaluateForce eval_vdw_force = EvaluateForce::NO,
                                  int system_index = 0, double clash_minimum_distance = 0.0,
                                  double clash_ratio = 0.0);

double2 evaluateAttenuated14Terms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                                  EvaluateForce eval_elec_force = EvaluateForce::NO,
                                  EvaluateForce eval_vdw_force = EvaluateForce::NO,
                                  int system_index = 0, double clash_minimum_distance = 0.0,
                                  double clash_ratio = 0.0);

double2 evaluateAttenuated14Terms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                                  EvaluateForce eval_elec_force = EvaluateForce::NO,
                                  EvaluateForce eval_vdw_force = EvaluateForce::NO,
                                  int system_index = 0, double clash_minimum_distance = 0.0,
                                  double clash_ratio = 0.0);

double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  const CoordinateFrameReader cfr, ScoreCard *ecard,
                                  int system_index = 0, double clash_minimum_distance = 0.0,
                                  double clash_ratio = 0.0);

double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  const CoordinateFrameWriter &cfw, ScoreCard *ecard,
                                  int system_index = 0, double clash_minimum_distance = 0.0,
                                  double clash_ratio = 0.0);

double2 evaluateAttenuated14Terms(const AtomGraph &ag, const CoordinateFrame &cf,
                                  ScoreCard *ecard, int system_index = 0,
                                  double clash_minimum_distance = 0.0, double clash_ratio = 0.0);

double2 evaluateAttenuated14Terms(const AtomGraph *ag, const CoordinateFrame &cf,
                                  ScoreCard *ecard, int system_index = 0,
                                  double clash_minimum_distance = 0.0, double clash_ratio = 0.0);

template <typename Tcoord, typename Tcalc>
double2 evaluateAttenuated14Terms(const ValenceKit<Tcalc> vk,
                                  const CoordinateSeriesReader<Tcoord> csr, ScoreCard *ecard,
                                  int system_index = 0, int force_scale_bits = 23,
                                  Tcalc clash_minimum_distance = 0.0, Tcalc clash_ratio = 0.0);

template <typename Tcoord, typename Tcalc>
double2 evaluateAttenuated14Terms(const ValenceKit<Tcalc> vk,
                                  const CoordinateSeriesWriter<Tcoord> csr, ScoreCard *ecard,
                                  int system_index = 0, int force_scale_bits = 23,
                                  Tcalc clash_minimum_distance = 0.0, Tcalc clash_ratio = 0.0);
/// \}

/// \brief Evaluate a positional restraint.  Return the energy penalty.
///
/// \param p_atom           The atom being restraint (in distance, angle, or dihedral restraints
///                         the atoms become I, J, K, and L)
/// \param step_number      The current step of the simulation, relevant only if the restraint
///                         evolves with time
/// \param init_step        Array of initial steps at which the restraint engages
/// \param finl_step        Array of final steps at which the restraint reaches its mature state
/// \param init_xy          Array of initial Cartesian X and Y target coordinates
/// \param finl_xy          Array of final Cartesian X and Y target coordinates
/// \param init_z           Array of initial Cartesian Z target coordinates
/// \param finl_z           Array of final Cartesian Z target coordinates
/// \param init_keq         Array of initial harmonic stiffnesses on the right and left sides
/// \param finl_keq         Array of final harmonic stiffnesses on the right and left sides
/// \param init_r           Array of initial displacement parameters
/// \param finl_r           Array of final displacement parameters
/// \param xcrd             Cartesian X coordinates of all particles
/// \param ycrd             Cartesian Y coordinates of all particles
/// \param zcrd             Cartesian Z coordinates of all particles
/// \param umat             Box space transformation matrix
/// \param invu             Inverse transformation, fractional coordinates back to real space
/// \param unit_cell        The unit cell type, i.e. triclinic
/// \param xfrc             Cartesian X forces acting on all particles
/// \param yfrc             Cartesian Y forces acting on all particles
/// \param zfrc             Cartesian Z forces acting on all particles
/// \param eval_force       Flag to have forces due to the restraint evaluated
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
Tcalc evalPosnRestraint(int p_atom, int step_number, int init_step, int finl_step,
                        const Tcalc2 init_xy, const Tcalc2 finl_xy, Tcalc init_z, Tcalc finl_z,
                        const Tcalc2 init_keq, const Tcalc2 finl_keq, const Tcalc4 init_r,
                        const Tcalc4 finl_r, const Tcoord* xcrd, const Tcoord* ycrd,
                        const Tcoord* zcrd, const double* umat, const double* invu,
                        UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
                        EvaluateForce eval_force, Tcalc inv_gpos_factor = 1.0,
                        Tcalc force_factor = 1.0);

/// \brief Evaluate distance restraints.  See the descriptions in evalPosnRestraint (above) for
///        explanations of each input parameter.  The restraint energy is returned.
///
/// Parameters for this function follow eponymous inputs to evalPosnRestraint(), above.
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
Tcalc evalBondRestraint(int i_atom, int j_atom, int step_number, int init_step, int finl_step,
                        const Tcalc2 init_keq, const Tcalc2 finl_keq, const Tcalc4 init_r,
                        const Tcalc4 finl_r, const Tcoord* xcrd, const Tcoord* ycrd,
                        const Tcoord* zcrd, const double* umat, const double* invu,
                        UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
                        EvaluateForce eval_force, Tcalc inv_gpos_factor = 1.0,
                        Tcalc force_factor = 1.0);

/// \brief Evaluate three-point angle restraints.  See the descriptions in evalPosnRestraint
///        (above) for explanations of each input parameter.  The restraint energy is returned.
///
/// Parameters for this function follow eponymous inputs to evalPosnRestraint(), above.
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
Tcalc evalAnglRestraint(int i_atom, int j_atom, int k_atom, int step_number, const int init_step,
                        const int finl_step, const Tcalc2 init_keq, const Tcalc2 finl_keq,
                        const Tcalc4 init_r, const Tcalc4 finl_r, const Tcoord* xcrd,
                        const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                        const double* invu, UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                        Tforce* zfrc, EvaluateForce eval_force, Tcalc inv_gpos_factor = 1.0,
                        Tcalc force_factor = 1.0);

/// \brief Evaluate four-point dihedral restraints.  See the descriptions in evalPosnRestraint
///        (above) for explanations of each input parameter.  The restraint energy is returned.
///
/// Parameters for this function follow eponymous inputs to evalPosnRestraint(), above.
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
Tcalc evalDiheRestraint(int i_atom, int j_atom, int k_atom, int l_atom, int step_number,
                        const int init_step, const int finl_step, const Tcalc2 init_keq,
                        const Tcalc2 finl_keq, const Tcalc4 init_r, const Tcalc4 finl_r,
                        const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                        const double* umat, const double* invu, UnitCellType unit_cell,
                        Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, EvaluateForce eval_force,
                        Tcalc inv_gpos_factor = 1.0, Tcalc force_factor = 1.0);

/// \brief Evaluate flat-bottom bimodal harmonic restraints of the form used in Amber's sander and
///        pmemd programs for NMR annealing calculations.  Time dependence of the restraints is
///        recognized, although there is not the same diversity of time evolution functions in
///        STORMM.
///
/// Overloaded:
///   - Evaluate based on a PhaseSpace object, with the option to compute and store forces
///   - Evaluate energy only based on a CoordinateFrame abstract
///   - Pass the restraint collection by pointer, by reference, or the reader abstract by value
///
/// \param ra            Restraint apparatus applicable to the topology that describes the system
/// \param rar           Double-precision reader for the restraint apparatus 
/// \param ps            Coordinates, box size, and force accumulators (modified by this function)
/// \param psw           Coordinates, box size, and force accumulators (modified by this function)
/// \param cfr           Coordinates of all particles, plus box dimensions (if needed)
/// \param cfw           Coordinates of all particles, plus box dimensions (if needed)
/// \param ecard         Energy components and other state variables (volume, temperature, etc.)
///                      (modified by this function)
/// \param eval_force    Flag to have forces also evaluated
/// \param system_index  Index of the system to which this energy contributes
/// \param step_number   The step number at which the energy is being evaluated (may determine the
///                      restraint parameters by mixing their endpoint values)
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
double evaluateRestraints(const RestraintKit<Tcalc, Tcalc2, Tcalc4> rar, const Tcoord* xcrd,
                          const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                          const double* invu, const UnitCellType unit_cell, Tforce* xfrc,
                          Tforce* yfrc, Tforce* zfrc, ScoreCard *ecard,
                          EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0,
                          int step_number = 0, Tcalc inv_gpos_factor = 1.0,
                          Tcalc force_factor = 1.0);

double evaluateRestraints(const RestraintKit<double, double2, double4> rar,
                          PhaseSpaceWriter psw, ScoreCard *ecard,
                          EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0,
                          int step_number = 0);

double evaluateRestraints(const RestraintApparatus &ra, PhaseSpace *ps, ScoreCard *ecard,
                          EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0,
                          int step_number = 0);

double evaluateRestraints(const RestraintApparatus *ra, PhaseSpace *ps, ScoreCard *ecard,
                          EvaluateForce eval_force = EvaluateForce::NO, int system_index = 0,
                          int step_number = 0);

double evaluateRestraints(const RestraintKit<double, double2, double4> rar,
                          const CoordinateFrameReader cfr, ScoreCard *ecard, int system_index = 0,
                          int step_number = 0);

double evaluateRestraints(const RestraintKit<double, double2, double4> rar,
                          const CoordinateFrameWriter &cfw, ScoreCard *ecard, int system_index = 0,
                          int step_number = 0);

double evaluateRestraints(const RestraintApparatus &ra, const CoordinateFrameReader cfr,
                          ScoreCard *ecard, int system_index = 0, int step_number = 0);

double evaluateRestraints(const RestraintApparatus *ra, const CoordinateFrameReader cfr,
                          ScoreCard *ecard, int system_index = 0, int step_number = 0);

template <typename Tcoord, typename Tcalc, typename Tcalc2, typename Tcalc4>
double evaluateRestraints(const RestraintKit<Tcalc, Tcalc2, Tcalc4> rar,
                          const CoordinateSeriesReader<Tcoord> csr, ScoreCard *ecard,
                          int system_index = 0, int step_number = 0, int force_scale_bits = 23);

template <typename Tcoord, typename Tcalc, typename Tcalc2, typename Tcalc4>
double evaluateRestraints(const RestraintKit<Tcalc, Tcalc2, Tcalc4> rar,
                          const CoordinateSeriesWriter<Tcoord> csr, ScoreCard *ecard,
                          int system_index = 0, int step_number = 0, int force_scale_bits = 23);
/// \}
  
} // namespace energy
} // namespace stormm

#include "valence_potential.tpp"

#endif
