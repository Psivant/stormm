// -*-c++-*-
#ifndef STORMM_TRAJECTORY_TRIM_H
#define STORMM_TRAJECTORY_TRIM_H

#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Constants/behavior.h"
#include "DataTypes/common_types.h"
#include "Math/formulas.h"
#include "Math/matrix_ops.h"
#include "Math/vector_ops.h"
#include "Numerics/split_fixed_precision.h"
#include "Parsing/parsing_enumerators.h"
#include "Parsing/parse.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_enumerators.h"
#include "motion_sweeper.h"
#include "phasespace.h"

namespace stormm {
namespace trajectory {

using card::GpuDetails;
using constants::ExceptionResponse;
using constants::PrecisionModel;
using data_types::isSignedIntegralScalarType;
using parse::NumberFormat;
using parse::realToString;
using stmath::angleOnAxes;
using stmath::angleOnAxesf;
using stmath::crossProduct;
using stmath::invertSquareMatrix;
using stmath::leibnizDeterminant;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisReader;
using synthesis::PsSynthesisWriter;
using synthesis::SyAtomUpdateKit;
using topology::AtomGraph;

/// \brief Check that the inertial tensor is invertible, returning TRUE if so but FALSE otherwise.
///        An error message or warning may be produced, depending on the specified behavior.
///
/// \param inrt    The inertial tensor to examine
/// \param natom   The number of atoms in the system (for error reporting purposes)
/// \param policy  Indicate the course of action if the moment of inertia cannot be computed
///                due to a singular matrix (a colinear arrangement of atoms will do this)
template <typename T>
bool checkInertialTensor(const std::vector<T> &inrt, int natom,
                         ExceptionResponse policy = ExceptionResponse::WARN);

/// \brief Remove net translational and, if applicable, angular momentum from the system.  This
///        function will also move the center of mass of a non-periodic system (no boundary
///        conditions) to the origin of the coordinate system.
///
/// \param xcrd        Cartesian X positions of all particles
/// \param ycrd        Cartesian Y positions of all particles
/// \param zcrd        Cartesian Z positions of all particles
/// \param xcrd_ovrf   Extended precision bits for Cartesian X positions of all particles
/// \param ycrd_ovrf   Extended precision bits for Cartesian Y positions of all particles
/// \param zcrd_ovrf   Extended precision bits for Cartesian Z positions of all particles
/// \param xvel        Cartesian X velocities of all particles
/// \param yvel        Cartesian Y velocities of all particles
/// \param zvel        Cartesian Z velocities of all particles
/// \param xvel_ovrf   Extended precision bits for Cartesian X velocities of all particles
/// \param yvel_ovrf   Extended precision bits for Cartesian Y velocities of all particles
/// \param zvel_ovrf   Extended precision bits for Cartesian Z velocities of all particles
/// \param masses      Masses of all particles, in atomic units (Daltons, g/mol)
/// \param natom       The number of particles
/// \param gpos_scale  Scaling factor for global positions, if expressed in fixed precision
/// \param vel_scale   Scaling factor for particle velocities, if expressed in fixed precision
/// \param prec        The precision model to use in computing momenta
/// \param ps
/// \param poly_ps
/// \param ag
/// \param poly_ag
/// \param mos
/// \param gpu         Details of the GPU that will carry out calculations
/// \param policy      Indicate the course of action if the moment of inertia cannot be computed
///                    due to a singular matrix (a colinear arrangement of atoms will do this)
/// \{
template <typename Tcoord, typename Tmass, typename Tcalc>
void removeMomentum(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, int* xcrd_ovrf, int* ycrd_ovrf,
                    int* zcrd_ovrf, Tcoord* xvel, Tcoord* yvel, Tcoord* zvel, int* xvel_ovrf,
                    int* yvel_ovrf, int* zvel_ovrf, const Tmass* masses, UnitCellType unit_cell,
                    int natom, Tcalc gpos_scale = 1.0, Tcalc vel_scale = 1.0,
                    ExceptionResponse policy = ExceptionResponse::WARN);

void removeMomentum(PhaseSpace *ps, const AtomGraph *ag,
                    ExceptionResponse policy = ExceptionResponse::WARN);

void removeMomentum(PhaseSpace *ps, const AtomGraph &ag,
                    ExceptionResponse policy = ExceptionResponse::WARN);

void removeMomentum(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                    const PrecisionModel prec, ExceptionResponse policy = ExceptionResponse::WARN);

void removeMomentum(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                    const PrecisionModel prec, ExceptionResponse policy = ExceptionResponse::WARN);

void removeMomentum(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                    MotionSweeper *mos, const GpuDetails &gpu = null_gpu,
                    ExceptionResponse policy = ExceptionResponse::WARN);
/// \}

/// \brief Accumulate the center of mass and translational momentum in a fixed-precision framework.
///        This sets the stage for removal of the center of mass's motion.
///
/// \param mosw      Holds accumulators for the center of mass and translational momentum of each
///                  system
/// \param poly_auk  One of the abstracts produced by the topology synthesis, containing atomic
///                  masses
/// \param poly_psr  Read-only abstract of the coordinate synthesis, containing coordinates and
///                  velocities
/// \param gpu       Specifications of the GPU that will perform the calculations
void accumulateCenterOfMassMotion(MotionSweepWriter *mosw,
                                  const SyAtomUpdateKit<double, double2, double4> &poly_auk,
                                  const PsSynthesisReader &poly_psr,
                                  const GpuDetails &gpu = null_gpu);

/// \brief Recenter all systems in the synthesis based on the accumulated centers of mass, and
///        remove the net translational velocity.
///
/// \param poly_psw  Writeable abstract for the coordinate synthesis, containing coordinates to
///                  modify in the arrays it presents as "current"
/// \param mosr      Read-only abstract of the motion sweeper, containing accumulated (but
///                  non-normalized) centers of mass
/// \param gpu       Specifications of the GPU that will perform the calculations
void removeCenterOfMassMotion(PsSynthesisWriter *poly_psw, const MotionSweepReader &mosr,
                              const GpuDetails &gpu = null_gpu);

/// \brief Accumulate the total rotational momentum and moment of inertia for each system in the
///        synthesis.  The moment of inertia and rotational momentum calculations assume that each
///        system has been positioned with its center of mass on the origin and net velocity
///        removed.  Descriptions of input parameters follow from the function
///        accumulateCenterOfMassMotion(), above.
void accumulateAngularMomentum(MotionSweepWriter *mosw,
                               const SyAtomUpdateKit<double, double2, double4> &poly_auk,
                               const PsSynthesisReader &poly_psr,
                               const GpuDetails &gpu = null_gpu);

/// \brief Subtract off the net translational as well as angular velocity (if appropriate for a
///        system with no boundary conditions), given the pre-computed momenta and inertial moment
///        tensors for all systems in a synthesis.  Descriptions of input parameters follow from
///        removeCenterOfMassMotion(), above, in addition to:
///
/// \param policy  Indicate the course of action if the moment of inertia cannot be computed due
///                to a singular matrix
void removeAngularMomentum(PsSynthesisWriter *poly_psw, const MotionSweepReader &mosr,
                           const GpuDetails &gpu = null_gpu,
                           ExceptionResponse policy = ExceptionResponse::WARN);
  
} // namespace trajectory
} // namespace stormm

#include "trim.tpp"

#endif
