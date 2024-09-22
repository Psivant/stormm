// -*-c++-*-
#ifndef STORMM_MAP_DENSITY_H
#define STORMM_MAP_DENSITY_H

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

/// \brief Check that the themes of the density-collecting particle-mesh interaction grid and the
///        atom-bearing cell grid are compatible.
///
/// \param pm_theme        Type of density to be mapped onto the particle-mesh interaction grid
/// \param cg_theme        Type of non-bonded property which the cell grid is set up to track
void matchThemes(NonbondedTheme pm_theme, NonbondedTheme cg_theme);
  
/// \brief Compute the grid position at which to being mapping density (the B-spline coefficients
///        will guide contributions along all three grid axes in reverse order).  The index and
///        fractional offsets of the particle along all three dimensions will be computed.
///
/// \param x               Cartesian X coordinate of the particle
/// \param y               Cartesian Y coordinate of the particle
/// \param z               Cartesian Z coordinate of the particle
/// \param inv_lpos_scale  Inverse scaling factor for the coordinates.  Providing 1.0 will indicate
///                        that the particle coordinates x, y, and z are given as real numbers.
/// \param cg_mesh_ticks   The number of particle-mesh interaction grid points per side of each
///                        spatial decomposition cell
/// \param cell_i          Spatial decomposition cell index along the system's A axis
/// \param cell_j          Spatial decomposition cell index along the system's B axis
/// \param cell_k          Spatial decomposition cell index along the system's C axis
/// \param a_cof           Fractional offset of the particle along the mesh A axis (returned)
/// \param b_cof           Fractional offset of the particle along the mesh B axis (returned)
/// \param c_cof           Fractional offset of the particle along the mesh C axis (returned)
/// \param bspline_order   The order of the B-splines to compute and the trusted length of a_cof,
///                        b_cof, and c_cof
/// \param grid_a          Grid point along the A axis at which to begin mapping density (returned)
/// \param grid_b          Grid point along the B axis at which to begin mapping density (returned)
/// \param grid_c          Grid point along the C axis at which to begin mapping density (returned)
template <typename Tcalc, typename Tgrid>
void particleAlignment(Tcalc x, Tcalc y, Tcalc z, Tcalc inv_lpos_scale, const Tcalc* umat,
                       int cg_mesh_ticks, int cell_i, int cell_j, int cell_k, Tgrid *a_cof,
                       Tgrid *b_cof, Tgrid *c_cof, int bspline_order, int *grid_a, int *grid_b,
                       int *grid_c);

/// \brief Spread the density of a particle to the grid based on starting indices along all three
///        axes plus pre-computed B-spline coefficients.  The calculation and accumulation types
///        for this function depend on the PMIGrid object, not the CellGrid object.
///
/// \param a_cof                B-spline coefficients for spreading the particle's density along
///                             the system's A axis
/// \param b_cof                B-spline coefficients for spreading the particle's density along
///                             the system's B axis
/// \param c_cof                B-spline coefficients for spreading the particle's density along
///                             the system's C axis
/// \param bspline_order        The order of B-spline interpolation to the mesh, and the trusted
///                             length of a_cof, b_cof, and c_cof
/// \param grid_root_a          Root grid element for spreading density along the system's A axis
/// \param grid_root_b          Root grid element for spreading density along the system's B axis
/// \param grid_root_c          Root grid element for spreading density along the system's C axis
/// \param grid_dims            Dimensions of the grids in each system, with the dimensions along
///                             the A, B, and C axes given in the "x", "y", and "z" members of each
///                             tuple.  The offset for reading / writing system j's grid is given
///                             in the "w" member of grid_dims[j].
/// \param grid_data            The primary data array for the grid density representation, or the
///                             only array if the grid keeps real-valued density.
/// \param overflow             Overflow bits for fixed-precision grid accumulation.  Leaving this
///                             as nullptr will indicate that the grid handles only real numbers
///                             and does not support fixed-precision accumulation.
/// \param mesh_scaling_factor  The scaling factor applied to density contributions in
///                             fixed-precision accumulation on the mesh
template <typename Tcalc, typename Tgrid>
void spreadDensity(const Tcalc* a_cof, const Tcalc* b_cof, const Tcalc* c_cof, int bspline_order,
                   int grid_root_a, int grid_root_b, int grid_root_c, const uint4 grid_dims,
                   FFTMode fft_staging, Tgrid* grid_data, int* overflow = nullptr,
                   Tcalc mesh_scaling_factor = 1.0);

/// \brief Accumulate the density for particles of a specific system in one of its spatial
///        decomposition cells.  This highly templated function serves a number of precision
///        models.
///
/// Overloaded:
///   - Operate on a set of real-valued particle-mesh interaction grids
///   - Operate on a set of fixed-precision particle-mesh interaction grids
///
/// \param pm_wrt  Writeable abstract for real-valued particle-mesh interaction grids
/// \param pm_acc  Writeable abstract for fixed-precision particle-mesh interaction grids
/// \param sysid   Index of the system of interest, in the synthesis as well as in the cell grid
/// \param cell_i  Index of the cell to operate on along the system grid's A axis
/// \param cell_j  Index of the cell to operate on along the system grid's B axis
/// \param cell_k  Index of the cell to operate on along the system grid's C axis
/// \param cgr     The cell grid, containing coordinates in local axes and the means to determine
///                their alignment to their respective system grids
/// \param synbk   Non-bonded parameter tables for the synthesis of all systems
/// \{
template <typename T, typename Tacc, typename Tcalc, typename Tcalc2, typename T4>
void accumulateCellDensity(PMIGridWriter *pm_wrt, int sysid, int cell_i, int cell_j, int cell_k,
                           const CellGridReader<T, Tacc, Tcalc, T4> &cgr,
                           const SyNonbondedKit<Tcalc, Tcalc2> &synbk);

template <typename T, typename Tacc, typename Tcalc, typename Tcalc2, typename T4>
void accumulateCellDensity(PMIGridAccumulator *pm_acc, int sysid, int cell_i, int cell_j,
                           int cell_k, const CellGridReader<T, Tacc, Tcalc, T4> &cgr,
                           const SyNonbondedKit<Tcalc, Tcalc2> &synbk);
/// \}

/// \brief Given a synthesis of coordinates in a CellGrid and associated particle-mesh interaction
///        grids, compute the density (due to charge or dispersion sources).  The calculation mode
///        of the cell grid will determine the calculation mode of the mapping.
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
/// \param pm        The particle-mesh interaction grids
/// \param pm_acc    Accumulator abstract for the particle-mesh interaction grids
/// \param cg        The cell grids with localized coordinates
/// \param v_cgr     Template-less abstract for the cell grid (used to cross the C++ : HPC
///                  boundary)
/// \param poly_ag   The topology synthesis containing charge or dispersion parameters
/// \param synbk     Non-bonded parameter abstract from the topology synthesis, taken at the
///                  desired precision for the calculations and able to convey this information
///                  across the C++ : HPC boundary
/// \param launcher  Compendium of launch parameters for core MM-related kernels, including the
///                  particle-mesh density mapping methods
/// \param lp        Launch parameters for the appropriate HPC kernel (an "abstract" of launcher)
/// \param approach  Indicate a particular HPC method to use in mapping particle density to the
///                  particle-mesh interaction grids
/// \{
template <typename T, typename Tacc, typename Tcalc, typename T4>
void mapDensity(PMIGrid *pm, const CellGrid<T, Tacc, Tcalc, T4> *cg,
                const AtomGraphSynthesis *poly_ag);

template <typename T, typename Tacc, typename Tcalc, typename T4>
void mapDensity(PMIGrid *pm, const CellGrid<T, Tacc, Tcalc, T4> &cg,
                const AtomGraphSynthesis &poly_ag);

#ifdef STORMM_USE_HPC
void mapDensity(PMIGridWriter *pm_wrt, PMIGridAccumulator *pm_acc, MMControlKit<double> *ctrl,
                const CellGridReader<void, void, void, void> &v_cgr, size_t cg_tmat,
                const SyNonbondedKit<double, double2> &synbk, int block_count, const int2 lp,
                QMapMethod approach, PMIGrid *pm);

void mapDensity(PMIGridWriter *pm_wrt, PMIGridAccumulator *pm_acc, MMControlKit<float> *ctrl,
                const CellGridReader<void, void, void, void> &v_cgr, size_t cg_tmat,
                const SyNonbondedKit<float, float2> &synbk, int block_count, const int2 lp,
                QMapMethod approach, PMIGrid *pm);

template <typename T, typename Tacc, typename Tcalc, typename T4>
void mapDensity(PMIGrid *pm, MolecularMechanicsControls *mm_ctrl,
                const CellGrid<T, Tacc, Tcalc, T4> *cg, const AtomGraphSynthesis *poly_ag,
                const CoreKlManager &launcher, QMapMethod approach);
#endif

void mapDensity(PMIGrid *pm, const AtomGraphSynthesis *poly_ag);

void mapDensity(PMIGrid *pm, const AtomGraphSynthesis &poly_ag);

template <typename Tcalc>
std::vector<Tcalc> mapDensity(const CoordinateFrameReader &cfr, const NonbondedKit<Tcalc> &nbk,
                              NonbondedTheme theme, FFTMode fft_staging, int grid_dim_a,
                              int grid_dim_b, int grid_dim_c, int order = default_bspline_order,
                              BSplineUnity unification = BSplineUnity::CENTER_FILL);

std::vector<double> mapDensity(const CoordinateFrame *cf, const AtomGraph *ag,
                               NonbondedTheme theme, FFTMode fft_staging, int grid_dim_a = -1,
                               int grid_dim_b = -1, int grid_dim_c = -1,
                               int order = default_bspline_order);

std::vector<double> mapDensity(const CoordinateFrame &cf, const AtomGraph &ag,
                               NonbondedTheme theme, FFTMode fft_staging, int grid_dim_a = -1,
                               int grid_dim_b = -1, int grid_dim_c = -1,
                               int order = default_bspline_order);
/// \}

/// \brief Unroll the call to an appropriately templated mapDensity function at the level of the
///        accumulator.  Unrolling at the level of cell dimension matrix / coordinate
///        representations is done in the non-templated mapDensity overloads.
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
void unrollMapDensityCall(PMIGrid *pm, size_t cg_tacc, size_t cg_tcalc,
                          const AtomGraphSynthesis *poly_ag);

template <typename T, typename Tacc, typename T4>
void unrollMapDensityCall(PMIGrid *pm, size_t cg_tcalc, const AtomGraphSynthesis *poly_ag);
/// \}

#ifdef STORMM_USE_HPC
/// \brief Launch the appropriate kernel to handle density mapping.  Like other template handoffs,
///        this will accept a template-free abstract from the C++ code and restore its
///        characteristics on the HPC side.  The data types of the non-bonded parameter kit will
///        be chosen based on the precision level of the calculation, and convey this information
///        implicitly.
///
/// Overloaded:
///   - Provide a writeable, real-valued abstract of the particle-mesh interaction grid (this will
///     launch the optimized kernel, and raise an error if that kernel is unavailable)
///   - Provide a writeable, fixed-precision abstract of the particle-mesh interaction grid (this
///     will launch the general-purpose kernel applying atomicAdd() for all atoms)
///   - Perform calculations in single- or double-precision 
///
/// \param pm_wrt    Real-valued abstract of the particle-mesh interaction grid
/// \param pm_acc    Abstract of the particle-mesh interaction grid prepared for fixed-precision
///                  accumulation
/// \param overflow  Indicate whether fixed-precision accumulation can occur without including a
///                  32-bit overflow accumulator
/// \param ctrl      Tracking object containing the time step and global counters directing
///                  asynchronous work unit progress
/// \param v_cgr     Read-only, template-free abstract for the cell grid
/// \param cg_tmat   Representation of the spatial decomposition cell dimensions in the cell grid
///                  (transformation matrices from fractional coordinates in each cell into
///                  Cartesian space).  Accepted values include int, llint, float, and double.
/// \param synbk     Non-bonded parameters for all particles.  This conveys the 
/// \param lp        Launch parameters for the mapping kernel, must correspond to the order
///                  presented by pm_wrt and the precision model conveyed by synbk.  As in other
///                  contexts, the block and thread count tuple produced by the kernel manager is
///                  equivalent to the abstracts produced by objects of the underlying C++ classes,
///                  the PMIGrid and the CellGrid.  This function is designed to take pre-assembled
///                  abstracts for the fastest possible execution.
/// \{
void launchShrAccDensityKernel(PMIGridWriter *pm_wrt, bool overflow, MMControlKit<double> *ctrl,
                               const CellGridReader<void, void, void, void> &v_cgr, size_t cg_tmat,
                               const SyNonbondedKit<double, double2> &synbk, const int2 lp);
  
void launchShrAccDensityKernel(PMIGridWriter *pm_wrt, bool overflow, MMControlKit<float> *ctrl,
                               const CellGridReader<void, void, void, void> &v_cgr, size_t cg_tmat,
                               const SyNonbondedKit<float, float2> &synbk, const int2 lp);

void launchGenPrpDensityKernel(PMIGridAccumulator *pm_acc,
                               const CellGridReader<void, void, void, void> &v_cgr, size_t cg_tmat,
                               const SyNonbondedKit<double, double2> &synbk, const int2 lp);
  
void launchGenPrpDensityKernel(PMIGridAccumulator *pm_acc,
                               const CellGridReader<void, void, void, void> &v_cgr, size_t cg_tmat,
                               const SyNonbondedKit<float, float2> &synbk, const int2 lp);
/// \}
#endif

} // namespace energy
} // namespace stormm

#include "map_density.tpp"

#endif
