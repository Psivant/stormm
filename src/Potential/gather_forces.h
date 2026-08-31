// -*-c++-*-
#ifndef STORMM_GATHER_FORCES_H
#define STORMM_GATHER_FORCES_H

#include "copyright.h"
#include "Accelerator/core_kernel_manager.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/bspline.h"
#include "Math/math_enumerators.h"
#include "Math/rounding.h"
#include "MolecularMechanics/mm_controls.h"
#include "Numerics/split_fixed_precision.h"
#include "Structure/local_arrangement.h"
#include "Structure/structure_enumerators.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "cellgrid.h"
#include "energy_enumerators.h"
#include "map_density.h"
#include "pmigrid.h"

namespace stormm {
namespace energy {

using card::CoreKlManager;
using data_types::isFloatingPointScalarType;
using mm::MolecularMechanicsControls;
using mm::MMControlKit;
using numerics::hostDoubleToInt63;
using numerics::hostDoubleToInt95;
using stmath::bSpline;
using stmath::bSplineNoUnity;
using stmath::BSplineUnity;
using stmath::roundUp;
using structure::imageCoordinates;
using structure::ImagingMethod;
using synthesis::AtomGraphSynthesis;
using synthesis::SyNonbondedKit;
using topology::AtomGraph;
using topology::NonbondedKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;

/// \brief Interpolate forces onto a particle based on a grid-based potential and the derivatives
///        of B-spline coefficients with which that particle's density was mapped onto the grid.
///        The calculation and grid-based potential types are both expected to be real-valued.
///
/// \param a_cof          Particle B-spline coefficients along the unit cell A axis, with the
///                       particle's inherent weight (e.g. charge) folded in
/// \param b_cof          Particle B-spline coefficients along the unit cell B axis
/// \param c_cof          Particle B-spline coefficients along the unit cell C axis
/// \param da_cof         Derivatives of B-spline coefficients along the unit cell A axis
/// \param db_cof         Derivatives of B-spline coefficients along the unit cell B axis
/// \param dc_cof         Derivatives of B-spline coefficients along the unit cell C axis
/// \param bspline_order  The order of B-spline interpolation
/// \param pme_potential  Real-valued array containing the grid-based potential
/// \param grid_root_a    Root grid element for gathering density along the system's A axis
/// \param grid_root_b    Root grid element for gathering density along the system's B axis
/// \param grid_root_c    Root grid element for gathering density along the system's C axis
/// \param grid_dims      Dimensions of the system's grid, with the point counts along the A, B,
///                       and C axes in the "x", "y", and "z" members of the tuple.  The starting
///                       element along the grid_data variable (see below) is given in the tuple's
///                       "w" member.
/// \param fft_staging    Indicate whether grid_data is configured to handle in-place FFTs, in
///                       which case the A dimension will be padded to a number of elements
///                       2 * ((grid_dims.x / 2) + 1), or if out-of-place FFTs are being performed,
///                       in which case no grid padding is in effect.
template <typename Tgrid, typename Tcalc, typename Tacc>
void pullPMEForces(const Tgrid* a_cof, const Tgrid* b_cof, const Tgrid* c_cof, const Tgrid* da_cof,
                   const Tgrid* db_cof, const Tgrid* dc_cof, int bspline_order,
                   const Tgrid* pme_potential, int grid_root_a, int grid_root_b, int grid_root_c,
                   const uint4 grid_dims, FFTMode fft_staging, const Tcalc* umat,
                   uint cg_img_index, Tacc* cg_xfrc, Tacc* cg_yfrc, Tacc* cg_zfrc, int* cg_xovrf,
                   int* cg_yovrf, int* cg_zovrf, Tacc* cg_net_frc, int* cg_net_ovrf,
                   Tcalc cg_frc_scl);

/// \brief Interpolate forces onto particles within one spatial decomposition cell of a specific
///        system.
///
/// \param pmigr   Read-only abstract of the particle-mesh interaction grid, containing the
///                grid-based potential for each system
/// \param sysid   Index of the system of interest, in the synthesis as well as in the CellGrid
/// \param cell_i  Index of the cell to operate on along the system grid's A axis
/// \param cell_j  Index of the cell to operate on along the system grid's B axis
/// \param cell_k  Index of the cell to operate on along the system grid's C axis
/// \param cgw     Writeable abstract of the CellGrid, accepting force contributions in
///                fixed-precision format.  Much of the templating for this function is based on
///                the configuration of the CellGrid.
/// \param synbk   Non-bonded parameters from the topology synthesis to which the CellGrid is tied
template <typename T, typename Tacc, typename Tcalc, typename Tcalc2, typename T4>
void gatherCellForces(CellGridWriter<T, Tacc, Tcalc, T4> *cgr, const PMIGridReader &pmigr,
                      const PsSynthesisBorders &pssb, int sysid, int cell_i, int cell_j,
                      int cell_k, const SyNonbondedKit<Tcalc, Tcalc2> &synbk);

/// \brief Given a synthesis of systems with particles arranged in a CellGrid, potential fields
///        mapped in a particle interaction grid, and perhaps an available GPU, call the
///        appropriate interpolation protocol to extract forces from the mesh-based potential.
///
/// Overloaded:
///   - Provide the cell grid as an explicit parameter
///   - Extract the cell grid coordinates from the particle-mesh interaction grid object
///   - Provide the cell grid and topology objects by const pointer or by const reference
///   - Provide a single topology, coordinate frame, and an indication of the non-bonded property
///   - Provide appropriate abstracts of the cell grid and particle-mesh interaction grid objects,
///     plus launch parameters.  In this case, a pointer to the original particle-mesh interaction
///     grids is still provided, to ensure that the resulting format of the data is reflected in
///     the original object despite having circumvented the need to create new abstracts with each
///     call.
///
/// \param pm        Particle-mesh interaction grids for all systems in the synthesis
/// \param pmigr     Abstract of the particle-mesh interaction grids
/// \param cg        The cell grids, with localized coordinates
/// \param v_cgw     Template-free but writeable abstract for the cell grid.  This is used to cross
///                  the C++ : HPC boundary.
/// \param poly_ag   The synthesis of topologies for all systems, containing non-bonded parameters
///                  essential for calculating the force on each particle
/// \param synbk     Non-bonded parameter abstract from the topology synthesis, taken at the
///                  desired precision for the calculations and able to convey this information
///                  across the C++ : HPC boundary
/// \param launcher  Compendium of launch parameters for core MM-related kernels, including the
///                  particle-mesh force interpolation methods
/// \param lp        Launch parameters for the appropriate HPC kernel (an "abstract" of launcher)
/// \{
template <typename T, typename Tacc, typename Tnb_calc, typename Tnb_calc2, typename T4>
void gatherForces(CellGridWriter<T, Tacc, Tnb_calc, T4> *cgw, const PMIGridReader &pmigr,
                  const PsSynthesisBorders &pssb,
                  const SyNonbondedKit<Tnb_calc, Tnb_calc2> &synbk);

template <typename T, typename Tacc, typename Tcalc, typename T4>
void gatherForces(CellGrid<T, Tacc, Tcalc, T4> *cg, const PMIGrid *pm,
                  const AtomGraphSynthesis *poly_ag);

template <typename T, typename Tacc, typename Tcalc, typename T4>
void gatherForces(CellGrid<T, Tacc, Tcalc, T4> *cg, const PMIGrid &pm, 
                  const AtomGraphSynthesis &poly_ag);

#ifdef STORMM_USE_HPC
void gatherForces(CellGridWriter<void, void, void, void> *v_cgr, const PMIGridReader &pmigr,
                  MMControlKit<double> *ctrl, size_t cg_tmat,
                  const SyNonbondedKit<double, double2> &synbk, int block_count, const int2 lp,
                  const PMIGrid *pm);

void gatherForces(CellGridWriter<void, void, void, void> *v_cgr, const PMIGridReader &pmigr,
                  MMControlKit<float> *ctrl, size_t cg_tmat,
                  const SyNonbondedKit<float, float2> &synbk, int block_count, const int2 lp,
                  const PMIGrid *pm);

template <typename T, typename Tacc, typename Tcalc, typename T4>
void gatherForces(CellGrid<T, Tacc, Tcalc, T4> *cg, const PMIGrid *pm,
                  MolecularMechanicsControls *mm_ctrl, const AtomGraphSynthesis *poly_ag,
                  const CoreKlManager &launcher);

template <typename T, typename Tacc, typename Tcalc, typename T4>
void gatherForces(CellGrid<T, Tacc, Tcalc, T4> *cg, const PMIGrid &pm,
                  MolecularMechanicsControls *mm_ctrl, const AtomGraphSynthesis &poly_ag,
                  const CoreKlManager &launcher);
#endif

void gatherForces(PMIGrid *pm, const AtomGraphSynthesis *poly_ag);

void gatherForces(PMIGrid *pm, const AtomGraphSynthesis &poly_ag);

template <typename Tdata, typename Tcalc, typename Tcalc3>
std::vector<Tcalc3> gatherForces(const CoordinateFrameReader &cfr, const Tdata* potential_field,
                                 const NonbondedKit<Tcalc> &nbk, NonbondedTheme theme,
                                 FFTMode fft_staging, int grid_dim_a, int grid_dim_b,
                                 int grid_dim_c, int order = default_bspline_order,
                                 BSplineUnity unification = BSplineUnity::CENTER_FILL);

template <typename Tdata, typename Tcalc>
std::vector<double3> gatherForces(const CoordinateFrame *cf, const Tdata* potential_field,
                                  const AtomGraph *ag, NonbondedTheme theme, FFTMode fft_staging,
                                  int grid_dim_a = -1, int grid_dim_b = -1, int grid_dim_c = -1,
                                  int order = default_bspline_order,
                                  PrecisionModel prec = PrecisionModel::DOUBLE,
                                  BSplineUnity unification = BSplineUnity::CENTER_FILL);

template <typename Tdata, typename Tcalc>
std::vector<double3> gatherForces(const CoordinateFrame &cf, const Tdata* potential_field,
                                  const AtomGraph &ag, NonbondedTheme theme, FFTMode fft_staging,
                                  int grid_dim_a = -1, int grid_dim_b = -1, int grid_dim_c = -1,
                                  int order = default_bspline_order,
                                  PrecisionModel prec = PrecisionModel::DOUBLE,
                                  BSplineUnity unification = BSplineUnity::CENTER_FILL);
/// \}

/// \brief Unroll the call to an appropriately templated gatherForces function, recovering some
///        template parameters of the cell grid object underlying a particle-mesh interaction grid.
///
/// Overloaded:
///   - Unroll the call to an appropriately templated overload of this function at the level of
///     the accumulator.
///   - Unroll the call to an appropriately templated mapDensity function at the level of the
///     the calculation mode.
///
/// \param pm        The particle-mesh interaction grids
/// \param cg_tacc   Detected data type ID for accumulation of forces in the cell grid
/// \param cg_tcalc  Detected data type ID for calculations on atoms in the cell grid
/// \param poly_ag   The topology synthesis containing charge or dispersion parameters
/// \{
template <typename T, typename T4>
void unrollGatherForcesCall(PMIGrid *pm, size_t cg_tacc, size_t cg_tcalc,
                            const AtomGraphSynthesis *poly_ag);

template <typename T, typename Tacc, typename T4>
void unrollGatherForcesCall(PMIGrid *pm, size_t cg_tcalc, const AtomGraphSynthesis *poly_ag);
/// \}

#ifdef STORMM_USE_HPC
/// \brief Launch the kernel to handle force interpolation from the transformed long-range
///        potential on the mesh.  This will accept a template-free abstract from the C++ layer,
///        then restore its data type characteristics on the HPC side.  As was the case with
///        density mapping kernel launch, the data types of the non-bonded parameter kit will be
///        chosen based on the precision level of the calculation.
///
/// Overloaded:
///   - Perform calculations in single- or double-precision
///
/// \param v_cgw    Writeable but template-less abstract of the cell grid neighbor list
/// \param pm_rdr   Read-only abstract of the particle-mesh interaction grid, containing the
///                 relevant non-bonded potential.
/// \param cg_tmat  Representation of the spatial decomposition cell dimensions in the cell grid
///                 (transformation matrices from fractional coordinates in each cell into
///                 Cartesian space).  Accepted values include int, llint, float, and double.
/// \param pssb     Contains unit cell sizes and transformation matrices, needed to properly scale
///                 forces on each particle
/// \param synbk    Non-bonded parameters for all particles.  This conveys the density property of
///                 each particle, whether charge or dispersion force strength, for every particle
///                 in each system of the synthesis covered by the cell grid and particle-mesh
///                 interaction grid.
/// \param lp       Launch parameters for the mapping kernel, corresponding to the order presented
///                 by pm_rdr and the precision model conveyed by synbk.  As in other contexts, the
///                 block and thread count tuple produced by the kernel manager is equivalent to
///                 the abstracts produced by objects of the underlying C++ classes, the PMIGrid
///                 and the CellGrid.  This function is designed to take pre-assembled abstracts
///                 for the fastest possible execution.
/// \{
void launchGenForceGatheringKernel(CellGridWriter<void, void, void, void> *v_cgw,
                                   const PMIGridReader &pm_rdr, const size_t cg_tmat,
                                   const PsSynthesisBorders &pssb,
                                   const SyNonbondedKit<double, double2> &synbk, const int2 lp);

void launchGenForceGatheringKernel(CellGridWriter<void, void, void, void> *v_cgw,
                                   const PMIGridReader &pm_rdr, const size_t cg_tmat,
                                   const PsSynthesisBorders &pssb,
                                   const SyNonbondedKit<float, float2> &synbk, const int2 lp);
/// \}
#endif
  
} // namespace energy
} // namespace stormm

#include "gather_forces.tpp"

#endif
