// -*-c++-*-
#ifndef STORMM_PMEGRID_H
#define STORMM_PMEGRID_H

#include <climits>
#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/math_enumerators.h"
#include "Numerics/split_fixed_precision.h"
#include "Synthesis/phasespace_synthesis.h"
#include "cellgrid.h"

namespace stormm {
namespace energy {

using card::GpuDetails;
using card::Hybrid;
using card::HybridKind;
using card::HybridTargetLevel;
using constants::CartesianDimension;
using constants::PrecisionModel;
using constants::UnitCellAxis;
using stmath::FFTMode;
using synthesis::PhaseSpaceSynthesis;

/// \brief Bounds on some settings for the particle-mesh interaction grid
/// \{
constexpr int minimum_qspread_fp_bits = 24;
constexpr int minimum_ljspread_fp_bits = 24;
constexpr int maximum_qspread_fp_bits_sp = 56;
constexpr int maximum_qspread_fp_bits_dp = 88;
constexpr int default_qspread_fp_bits_sp = 32;
constexpr int default_qspread_fp_bits_dp = 64;
constexpr int maximum_ljspread_fp_bits_sp = 50;
constexpr int maximum_ljspread_fp_bits_dp = 82;
constexpr int default_ljspread_fp_bits_sp = 24;
constexpr int default_ljspread_fp_bits_dp = 56;
constexpr int mapping_nonoverflow_bits_dp = 61;
constexpr int mapping_nonoverflow_bits_sp = 29;
constexpr int default_bspline_order = 5;
constexpr int density_mapping_wu_size = 32;
/// \}

/// \brief Default limits on the dimensions of the atom-bearing region in particle-mesh
///        interpolation work units.  The format of each name is "max" + [interpolation order] +
///        [s = SINGLE precision mode] + "_atom_bearing_" +
///        [region_adim = length, in spatial decomposition cells, along the unit cell A axis,
///         cross_section = area, in spatial decomposition cells, of the cross section along the
///                         B and C axes].  These numbers are set to keep within a limit of
///        48kB __shared__ memory per thread block for register-based accumulation and 24kB
///        __shared__ memory per thread block for __shared__-based accumulation.
/// \{
#ifdef STORMM_USE_CUDA
constexpr int max_shared_acc_atom_bearing_region_adim = 28;
#else
constexpr int max_shared_acc_atom_bearing_region_adim = 9;
#endif
constexpr int max4s_grid_mapping_volume = 48;
constexpr int max4s_atom_bearing_region_adim = 7;
constexpr int max4s_atom_bearing_cross_section = 16;
constexpr int max5s_grid_mapping_volume = 48;
constexpr int max5s_atom_bearing_region_adim = 7;
constexpr int max5s_atom_bearing_cross_section = 16;
constexpr int max6s_grid_mapping_volume = 48;
constexpr int max6s_atom_bearing_region_adim = 4;
constexpr int max6s_atom_bearing_cross_section = 12;
/// \}

/// \brief A writeable abstract for the Particle-Mesh Interaction Grid class.
struct PMIGridWriter {

  /// \brief The constructor take a straight list of inputs for all member variables.
  PMIGridWriter(NonbondedTheme theme_in, PrecisionModel mode_in, FFTMode fftm_in, int fp_bits_in,
                int nsys_in, int order_in, int wu_count_in, int max_grid_points_in,
                const uint4* dims_in, double* ddata_in, float* fdata_in,
                const uint* work_units_in);

  /// \brief As with other abstracts, the presence of one or more const members forbids definition
  ///        of the copy and move assignment operators, but with no pointers to repair the default
  ///        copy and move constructors apply.
  /// \{
  PMIGridWriter(const PMIGridWriter &original) = default;
  PMIGridWriter(PMIGridWriter &&original) = default;
  /// \}

  const NonbondedTheme theme;  ///< The non-bonded property mapped to the particle-mesh interaction
                               ///<   grid
  const PrecisionModel mode;   ///< The mode in which the object is to operate
  const FFTMode fftm;          ///< The layout of the FFT grid, implying a mode in which the FFT
                               ///<   will be carried out
  const float shacc_fp_scale;  ///< Scaling factor used in fixed-precision accumulation when
                               ///<   accumulating density in __shared__ memory.
  const int nsys;              ///< The number of systems in the synthesis served by this object,
                               ///<   and a trusted length for the dims array below
  const int order;             ///< The interpolation order for B-pline mapping of density to the
                               ///<   particle-mesh interaction grid
  const int wu_count;          ///< The number of work units built to guide special-purpose GPU
                               ///<   kernels.  When multiplied by 32, this is the trusted length
                               ///<   of work_units.
  const int max_grid_points;   ///< The naximum number of grid points that any one work unit deals
                               ///<   with
  const uint4* dims;           ///< Dimensions of each PMI grid.  A-, B-, and C-axis lengths are
                               ///<   found in the "x", "y", and "z" members of the tuple, with the
                               ///<   offset for the start of each grid in the "w" member.
  double* ddata;               ///< Double-precision real data, bounded by indices found in the
                               ///<   "w" member of dims
  float* fdata;                ///< Single-precision real data, bounded by indices found in the
                               ///<   "w" member of dims
  const uint* work_units;      ///< Array of 32-element work units to guide certain types of
                               ///<   accumulation kernels
};

/// \brief A read-only abstract for the Particle-Mesh Interaction Grid class.
struct PMIGridReader {

  /// \brief The constructor take a straight list of inputs for all member variables.
  PMIGridReader(NonbondedTheme theme_in, PrecisionModel mode_in, FFTMode fftm_in, int fp_bits_in,
                int nsys_in, int order_in, const uint4* dims_in, const double* ddata_in,
                const float* fdata_in);

  /// \brief As with some other readers, the writer can be used as a constructor argument.
  ///
  /// \param w  A writeable abstract for a pre-existing PMIGrid object
  /// \{
  PMIGridReader(const PMIGridWriter &w);
  PMIGridReader(const PMIGridWriter *w);
  /// \}
  
  /// \brief As with other abstracts, the presence of one or more const members forbids definition
  ///        of the copy and move assignment operators, but with no pointers to repair the default
  ///        copy and move constructors apply.
  /// \{
  PMIGridReader(const PMIGridReader &original) = default;
  PMIGridReader(PMIGridReader &&original) = default;
  /// \}

  const NonbondedTheme theme;  ///< The non-bonded property mapped to the particle-mesh interaction
                               ///<   grid
  const PrecisionModel mode;   ///< The mode in which the object is to operate
  const FFTMode fftm;          ///< The layout of the FFT grid, implying a mode in which the FFT
                               ///<   will be carried out
  const float shacc_fp_scale;  ///< Scaling factor used in fixed-precision accumulation when
                               ///<   accumulating density in __shared__ memory.
  const int nsys;              ///< The number of systems in the synthesis served by this object,
                               ///<   and a trusted length for the dims array below
  const int order;             ///< The interpolation order for B-pline mapping of density to the
                               ///<   particle-mesh interaction grid
  const uint4* dims;           ///< Dimensions of each PMI grid.  A-, B-, and C-axis lengths are
                               ///<   found in the "x", "y", and "z" members of the tuple, with the
                               ///<   offset for the start of each grid in the "w" member.
  const double* ddata;         ///< Double-precision real data, bounded by indices found in the
                               ///<   "w" member of dims
  const float* fdata;          ///< Single-precision real data, bounded by indices found in the
                               ///<   "w" member of dims
};
  
/// \brief A writeable abstract which reinterprets some pointers to enable split fixed-precision
///        accumulation within the PMIGrid object.
struct PMIGridAccumulator {

  /// \brief The constructor takes inputs from the object's member variables and re-interprets
  ///        some of them as necessary.
  PMIGridAccumulator(NonbondedTheme theme_in, PrecisionModel mode_in, FFTMode fftm_in,
                     bool use_overflow_in, int fp_bits_in, int nsys_in, int order_in,
                     int wu_count_in, const uint4* dims_in, double* ddata_in, float* fdata_in,
                     int* overflow_in, const uint* work_units_in);

  /// \brief As with other abstracts, the presence of one or more const members forbids definition
  ///        of the copy and move assignment operators, but with no pointers to repair the default
  ///        copy and move constructors apply.
  /// \{
  PMIGridAccumulator(const PMIGridAccumulator &original) = default;
  PMIGridAccumulator(PMIGridAccumulator &&original) = default;
  /// \}
  
  const NonbondedTheme theme;  ///< The non-bonded property mapped to the particle-mesh interaction
                               ///<   grid
  const PrecisionModel mode;   ///< The mode in which the object is to operate
  const FFTMode fftm;          ///< The layout of the FFT grid, implying a mode in which the FFT
                               ///<   will be carried out
  const bool use_overflow;     ///< Flag to indicate that overflow accumulators must be used
  const int fp_bits;           ///< The number of bits used in fixed-precision accumulation of
                               ///<   whatever density the grids describe.
  const float fp_scale;        ///< Scaling factor applied to all contributions for fixed-precision
                               ///<   accumulation
  const int nsys;              ///< The number of systems in the synthesis served by this object,
                               ///<   and a trusted length for the dims array below
  const int order;             ///< The interpolation order for B-pline mapping of density to the
                               ///<   particle-mesh interaction grid
  const int order_squared;     ///< Square of the interpolation order, for tracking progress in
                               ///<   naive mapping methods when multiple warps may work on the
                               ///<   same atom
  const int order_cubed;       ///< Cube of the interpolation order, for tracking progress in 
                               ///<   naive mapping methods when multiple warps may work on the
                               ///<   same atom
  const int wu_count;          ///< The number of work units built to guide special-purpose GPU
                               ///<   kernels.  When multiplied by 32, this is the trusted length
                               ///<   of work_units.
  const uint4* dims;           ///< Dimensions of each PMI grid.  A-, B-, and C-axis lengths are
                               ///<   found in the "x", "y", and "z" members of the tuple, with the
                               ///<   offset for the start of each grid in the "w" member.
  llint* lldata;               ///< Primary accumulator for 95-bit split fixed precision data
  int* idata;                  ///< Primary accumulator for 63-bit split fixed precision data
  int* overflow;               ///< Overflow data to complement either lldata or idata, creating
                               ///<   95-bit and 63-bit split fixed-precision accumulators,
                               ///<   respectively.
  const uint* work_units;      ///< Array of 32-element work units to guide certain types of
                               ///<   accumulation kernels
};

/// \brief A read-only abstract which can interpret the PMIGrid data in split fixed-precision
///        format.
struct PMIGridFPReader {

  /// \brief The constructor takes inputs from the object's member variables and re-interprets
  ///        some of them as necessary.
  PMIGridFPReader(NonbondedTheme theme_in, PrecisionModel mode_in, FFTMode fftm_in,
                  bool use_overflow_in, int fp_bits_in, int nsys_in, int order_in,
                  const uint4* dims_in, const double* ddata_in, const float* fdata_in,
                  const int* overflow_in);

  /// \brief As with some other readers, the writer can be used as a constructor argument.
  ///
  /// \param w  A fixed-precision accumulation abstract for a pre-existing PMIGrid object
  /// \{
  PMIGridFPReader(const PMIGridAccumulator &w);
  PMIGridFPReader(const PMIGridAccumulator *w);
  /// \}

  /// \brief As with other abstracts, the presence of one or more const members forbids definition
  ///        of the copy and move assignment operators, but with no pointers to repair the default
  ///        copy and move constructors apply.
  /// \{
  PMIGridFPReader(const PMIGridFPReader &original) = default;
  PMIGridFPReader(PMIGridFPReader &&original) = default;
  /// \}

  const NonbondedTheme theme;  ///< The non-bonded property mapped to the particle-mesh interaction
                               ///<   grid
  const PrecisionModel mode;   ///< The mode in which the object is to operate
  const FFTMode fftm;          ///< The layout of the FFT grid, implying a mode in which the FFT
                               ///<   will be carried out
  const bool use_overflow;     ///< Flag to indicate that overflow accumulators are relevant
  const int fp_bits;           ///< The number of bits used in fixed-precision accumulation of
                               ///<   whatever density the grids describe.
  const float fp_scale;        ///< Scaling factor applied to all contributions for fixed-precision
                               ///<   accumulation
  const int nsys;              ///< The number of systems in the synthesis served by this object,
                               ///<   and a trusted length for the dims array below
  const int order;             ///< The interpolation order for B-pline mapping of density to the
                               ///<   particle-mesh interaction grid
  const uint4* dims;           ///< Dimensions of each PMI grid.  A-, B-, and C-axis lengths are
                               ///<   found in the "x", "y", and "z" members of the tuple, with the
                               ///<   offset for the start of each grid in the "w" member.
  const llint* lldata;         ///< Primary accumulator for 95-bit split fixed precision data
  const int* idata;            ///< Primary accumulator for 63-bit split fixed precision data
  const int* overflow;         ///< Overflow data to complement either lldata or idata, creating
                               ///<   95-bit and 63-bit split fixed-precision accumulators,
                               ///<   respectively.
};

/// \brief An object to hold a series of meshes for accumulating density from condensed-phase
///        molecular systems and transforming it into potential, whether for charges and
///        electrostatics or dispersion sources and van-der Waals interactions.  This supports the
///        particle-mesh interactions (PMI) of the particle-particle / particle-mesh potential
///        partitioning common in several simulation techniques.
class PMIGrid {
public:

  /// \brief The constructor builds the object off of a CellGrid object, which is in turn linked
  ///        to a particular coordinate synthesis (PhaseSpaceSynthesis).
  /// \{
  template <typename T, typename Tacc, typename Tcalc, typename T4>
  PMIGrid(const CellGrid<T, Tacc, Tcalc, T4> *cg_in, NonbondedTheme theme_in,
          int b_spline_order_in = default_bspline_order,
          PrecisionModel mode_in = PrecisionModel::SINGLE,
          FFTMode fft_staging_in = FFTMode::IN_PLACE, int fp_accumulation_bits_in = 0,
          int shared_fp_accumulation_bits_in = -1, const GpuDetails &gpu = null_gpu,
          const QMapMethod work_unit_configuration_in = QMapMethod::GENERAL_PURPOSE);

  template <typename T, typename Tacc, typename Tcalc, typename T4>
  PMIGrid(const CellGrid<T, Tacc, Tcalc, T4> &cg_in, NonbondedTheme theme_in,
          int b_spline_order_in = default_bspline_order,
          PrecisionModel mode_in = PrecisionModel::SINGLE,
          FFTMode fft_staging_in = FFTMode::IN_PLACE, int fp_accumulation_bits_in = 0,
          int shared_fp_accumulation_bits_in = -1, const GpuDetails &gpu = null_gpu,
          const QMapMethod work_unit_configuration_in = QMapMethod::GENERAL_PURPOSE);
  /// \}

  /// \brief The choice to re-cast arrays in the abstracts rather than have special POINTER-kind
  ///        Hybrid objects in the PMIGrid itself allows the default copy and move constructors,
  ///        as well as the copy and move assignment operators, to remain valid.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object placed on the right hand side of the assignment statement
  /// \{
  PMIGrid(const PMIGrid &original) = default;
  PMIGrid(PMIGrid &&original) = default;
  PMIGrid& operator=(const PMIGrid &original) = default;
  PMIGrid& operator=(PMIGrid &&original) = default;
  /// \}

  /// \brief Get the non-bonded propery mapped onto the particle-mesh interaction grid.
  NonbondedTheme getTheme() const;
  
  /// \brief Get the mode in which the object is set to operate.
  PrecisionModel getMode() const;

  /// \brief Get the method for performing real-to-complex and complex-to-real FFTs.  This implies
  ///        how the grid data itself may be padded.
  FFTMode getFFTStaging() const;

  /// \brief Indicate whether fixed-precision accumulation is enabled.  This will test whether
  ///        fp_accumulation_bits is greater than zero.
  bool fixedPrecisionEnabled() const;

  /// \brief Get the recommended density mapping method for the current set of grids.
  QMapMethod getRecommendedMappingMethod() const;
  
  /// \brief Get the work unit configuration of the object.
  QMapMethod getWorkUnitConfiguration() const;
  
  /// \brief Indicate whether the data is present as real numbers or fixed-precision accumulators.
  ///        Even with fixed-precision enabled, this object may hold data as real numbers for
  ///        non-accumulating tasks.
  bool dataIsReal() const;
  
  /// \brief Get the number of bits used in scaling numbers for fixed precision accumulation.
  int getFixedPrecisionBits() const;

  /// \brief Get the number of bits used in scaling numbers for fixed precision accumulation in
  ///        __shared__ accumulation density mapping kernels.
  int getSharedFixedPrecisionBits() const;
  
  /// \brief Get the number of systems with grid managed by this object.
  int getSystemCount() const;

  /// \brief Get the interpolation order to be used in mapping particle density to the
  ///        particle-mesh interaction grid.
  int getInterpolationOrder() const;

  /// \brief Get the total capacity allocated to store grid data.
  size_t getTotalCapacity() const;
  
  /// \brief Obtain the dimensions of one of the grids in the object
  ///
  /// Overloaded:
  ///   - Get all dimensions, including the offset in the data arrays, of one system's PMI grid
  ///   - Get the dimension of one system's grid along a particular axis
  ///
  /// \param system_index  Index of the system of interest within the associated synthesis
  /// \param uc_axis       The unit cell axis of interest
  /// \{
  uint4 getGridDimensions(int system_index) const;
  int getGridDimensions(int system_index, UnitCellAxis uc_axis) const;
  int getGridDimensions(int system_index, CartesianDimension uc_axis) const;
  /// \}

  /// \brief Get the number of work units used by the optimized kernel.
  int getWorkUnitCount() const;

  /// \brief Report the number of grid points present in the largest work unit.
  int getLargestWorkUnitGridPoints() const;

  /// \brief Report whether the PMIGrid will use short format accumulation.
  bool shortFormatAccumulation() const;

  /// \brief Report whether the PMIGrid will need to use overflow accumulators.
  bool useOverflowAccumulation() const;

  /// \brief Get the abstract for this object containing its precision mode, grid dimensions, and
  ///        pointers to relevant data arrays.
  ///
  /// Overloaded:
  ///   - Get a writeable abstract from a mutable object.
  ///   - Get a read-only abstract from a const object.
  ///
  /// \param tier  Specify whether to obtain pointers on the CPU host or GPU device
  /// \{
  PMIGridWriter data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  const PMIGridReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Get the fixed-precision abstract for this object, re-interpreting either main data
  ///        array and adding an integer array of overflow accumulators for accumulation protected
  ///        against the non-associativity that breaks bitwise reproducibility for floating point
  ///        parallel computations.
  ///
  /// Overloaded:
  ///   - Get a writeable abstract from a mutable object.
  ///   - Get a read-only abstract from a const object.
  ///
  /// \param tier  Specify whether to obtain pointers on the CPU host or GPU device
  /// \{
  PMIGridAccumulator fpData(HybridTargetLevel tier = HybridTargetLevel::HOST);
  const PMIGridFPReader fpData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}
  
  /// \brief Get the pointer to the attached CellGrid object, which provides the basis for the
  ///        grids / meshes contained in this object.  This will restore the type of the object
  ///        but must be invoked in a manner consistent with its original typing.
  template <typename T, typename Tacc, typename Tcalc, typename T4>
  const CellGrid<T, Tacc, Tcalc, T4>* getCellGridPointer() const;

  /// \brief Get an abstract of the attached cell grid object, with all templating cast away.  This
  ///        exists because the templated types of the pointer stored in the PMIGrid object are
  ///        placeholders for its true types.  Out of an abundance of caution, a CellGrid pointer
  ///        with the proper data types for all four components is created so that the CellGrid's
  ///        public member function templateFreeData() can be invoked.
  ///
  /// \param tier  Specify whether to obtain pointers, void-casted or otherwise, on the CPU host or
  ///              GPU device
  const CellGridReader<void, void, void, void>
  getTemplateFreeCellGridReader(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  
  /// \brief Get the type ID of the attached CellGrid object's cell dimension matrix expressions.
  ///        This is the "typename T" seen in many CellGrid class function declarations, and
  ///        the root data type of the four-tuples in which its coordinates are expressed.
  ///        Possible types include float, double, int, and llint.
  size_t getCellGridMatrixTypeID() const;

  /// \brief Get the type ID of the attached CellGrid object's accumulators.  This is "typename
  ///        Tacc" in various CellGrid class function declarations and invocations.  Acceptable
  ///        types include int and llint.
  size_t getCellGridAccumulatorTypeID() const;

  /// \brief Get the type ID of the attached CellGrid object's fractional space coordinate
  ///        transformations, which sets a precedent for any calculations that will be done with
  ///        an object of the class.  This is "typename Tcalc" in various CellGrid class function
  ///        declarations and invocations.
  size_t getCellGridCalculationTypeID() const;

  /// \brief Get the type ID of the attached CellGrid object's coordinate data, a four-tuple of
  ///        its cell dimension matrix data.  Acceptable types include float, double, int, and
  ///        llint.  This is "typename T4" in most CellGrid-associated function declarations and
  ///        invocations, sometimes "typename Tcrd".
  size_t getCellGridCoordinateTypeID() const;
  
  /// \brief Get the total substance, particle density, mapped to all grid elements for a
  ///        particular system.
  ///
  /// \param system_index  The system of interest
  double getTotalOnGrid(int system_index) const;

  /// \brief Export the particle-mesh interaction grid asssociated with one system in the
  ///        synthesis.  The resulting linearized data reflects an ordering of the grid that is
  ///        identical to the way it appears in the PMIGrid itself.
  ///
  /// \param system_index  The system of interest
  std::vector<double> getGrid(int system_index) const;

  /// \brief Get a pointer to the underlying coordinate synthesis.
  const PhaseSpaceSynthesis* getCoordinateSynthesisPointer() const;
  
  /// \brief Get a pointer to the PMIGrid object itself.
  const PMIGrid* getSelfPointer() const;
  
#if STORMM_USE_HPC
  /// \brief Upload all data from the CPU host to the GPU device.
  void upload();

  /// \brief Download all data from the GPU device to the CPU host memory.
  void download();
#endif

  /// \brief Change the mode in which the object is set to operate.  This will leave bounds arrays
  ///        unchanged but blow away data as the object is reconfigured.
  ///
  /// \param mode_in  The chosen mode for the object
  void setMode(PrecisionModel mode_in);

  /// \brief Change the character of the data format.  If called with TRUE the object will report
  ///        holding real-valued data.
  ///
  /// \param real_in  Specify whether the object holds real (TRUE) or fixed-precision data (FALSE)
  void setRealDataFormat(bool real_in = true);
  
  /// \brief Set the object to enable or disable fixed-precision accumulation.
  ///
  /// \param fp_accumulation_bits_in  The number of bits used in fixed-precision accumulation,
  ///                                 indicating the power of two by which to scale numbers going
  ///                                 int each sum
  /// \param shared_fp_bits_in        The number of bits used in fixed-precision accumulation with
  ///                                 a separate space in __shared__ memory.  While this will be
  ///                                 set to be consistent with fp_accumulation_bits_in if
  ///                                 possible, there is a critical difference: this 
  ///                                 accumulation is done in the interest of going directly to
  void prepareFixedPrecisionModel(int fp_accumulation_bits_in, int shared_fp_bits_in);

  /// \brief Set the recommended mapping method based on the density of the grid and available
  ///        resources.
  ///
  /// \param gpu  Specifications of the GPU that will perform the work
  void setRecommendedMappingMethod(const GpuDetails &gpu);

  /// \brief Prepare a set of work units for the synthesis of systems.  These work units will be
  ///        tailored to the synthesis for which the topology was created, but can be changed given
  ///        a new GPU.
  ///
  /// Overloaded:
  ///   - Tailor work units to a chosen a mapping method
  ///   - Tailor work units to the recommended mapping method
  ///
  /// \param approach  The chosen work unit mapping method (if the grid density or parameters do
  ///                  do not comport with the method, an alternative will be selected)
  /// \param gpu       Details of the available GPU
  /// \{
  void prepareWorkUnits(QMapMethod approach, const GpuDetails &gpu);
  void prepareWorkUnits(const GpuDetails &gpu);
  /// \}
  
  /// \brief Initialize the particle-mesh interaction grids, respecting the object's precision and
  ///        accumulation modes.  If fixed-precision is enabled, this will set the data content to
  ///        register as no longer "real numbers" (not imaginary numbers, but integers).
  ///
  /// \param tier  Indicate whether to initialize the density grids at the level of the CPU host or
  ///              the GPU device
  /// \param gpu   Details of the available GPU, if the initialization focuses on device memory
  void initialize(HybridTargetLevel tier = HybridTargetLevel::HOST,
                  const GpuDetails &gpu = null_gpu);

  /// \brief Convert fixed-precision integer data in the object to real data.  There is no function
  ///        to convert real back to fixed-precision integer data: the process is to accumulate in
  ///        the one mode and then convert to the other when finished.
  ///
  /// \param tier  Indicate whether to initialize the density grids at the level of the CPU host or
  ///              the GPU device
  /// \param gpu   Details of the available GPU, if the initialization focuses on device memory
  void convertToReal(HybridTargetLevel tier = HybridTargetLevel::HOST,
                     const GpuDetails &gpu = null_gpu);
  
private:

  /// The non-bonded property mapped to the grid
  NonbondedTheme theme;
  
  /// Indicate whether the object will run Fast Fouriet Transforms (FFTs) in single- or
  /// double-precision arithmetic.
  PrecisionModel mode; 

  /// Indicate whether the object will perform real-to-complex and complex-to-real FFTs in-place
  /// or out-of-place.
  FFTMode fft_staging;
  
  /// The number of bits used in scaling for fixed-precision accumulation.  If set to zero,
  /// fixed-precision accumulation is disabled in the object, the overflow array is not allocated
  /// (or de-allocated if necessary), and the fp_(d,f)grid_stack pointers will be unset.
  int fp_accumulation_bits;

  /// The number of bits after the decimal to be used in __shared__ memory accumulation based on
  /// atomic contributions with fixed precision.  This will be set based on fp_accumulation_bits
  /// if that parameter is set to anything other than zero, but can be set independently if that
  /// parameter is zero.
  int shared_fp_accumulation_bits;
  
  /// Flag to indicate whether the object currently holds real-value or fixed-precision data.
  /// This is initialized when the fixed-precision bit count is set, and toggled by calls to
  /// initialize() (for the beginning of accumulation) or convertToReal() (signifying the end of
  /// fixed-precision accumulation).
  bool data_is_real;

  /// Flag to indicate that the topology synthesis and fixed precision bit count used in any
  /// fixed-precision accumulation (whether in __shared__ memory or in __global__ memory) is
  /// compatible with methods that send results only to the primary accumulator.  A TRUE value
  /// in this flag will select kernels and functions that ignore possible overflow as this has
  /// been declared safe.
  bool use_short_format_accumulation;
  
  /// The number of systems with grids in the object.  This is taken from the coordinate synthesis
  /// responsible for creating the attached CellGrid.
  int system_count;

  /// The order of B-spline interpolation to be used when mapping particle density to the
  /// particle-mesh interaction grids
  int b_spline_order;

  /// The available __shared__ memory buffer size, in bytes, determined by examining the available
  /// GPU's specifications.
  int shared_acc_buffer_size;
  
  /// A recommended mapping method can be determined given a set of GPU resources.
  QMapMethod recommendation;

  /// Regardless of the recommended mapping method, the work units can be configured for a
  /// particular approach.
  QMapMethod work_unit_configuration;
  
  /// List of the grid dimensions for each system.  The dimensions along the A, B, and C axes are
  /// given in the "x", "y", and "z" members of the tuples.  The "w" member holds the overall
  /// offset of the mesh within the master data array.
  Hybrid<uint4> grid_dimensions;

  /// Length of the data array and (if applicable) the overflow data for fixed-precision
  /// accumulation.  The size must not exceed 2^32, or roughly 16GB of grid data amongst all
  /// systems, per the limits of indexing in grid_dimensions.
  size_t capacity;

  /// The number of work units to execute in the optimized kernel.
  int work_unit_count;

  /// The largest number of grid points involved in any one work unit
  int largest_work_unit_grid_points;
  
  /// The master array of FFT-ready double-precision data.  This can also serve as the primary
  /// accumulator for 95-bit fixed-precision accumulations that will be transformed into the
  /// real-valued data.  Allocated only if mode is set to DOUBLE.
  Hybrid<double> dgrid_stack;

  /// The master array of FFT-ready single-precision data.  This can also serve as the primary
  /// accumulator for 63-bit fixed-precision accumulations that will be transformed into the
  /// real-valued data.  Allocated only if mode is set to SINGLE.
  Hybrid<float> fgrid_stack;

  /// An auxiliary array of accumulators, the overflow bits, for fixed-precision work prior to
  /// converting data into real number format.  For either mode of operation, the primary data
  /// arrays will be the memory which, in another context, might be used for FFTs.
  Hybrid<int> overflow_stack;

  /// An array holding bitmask-based work units for filling out the particle-mesh interaction grids
  /// based on the cell grids.  Each work unit is a series of 32 numbers (independent of the
  /// architecture, although this is NVIDIA's warp size as of CUDA 12 and has always been).  The
  /// numbers in each work unit must be interpreted based on the work_unit_configuration setting
  /// (see above).  For ACC_SHARED, the work unit indices read:
  ///
  /// Index   Description
  ///   0     Index of the system to which the work unit pertains
  ///  1- 3   Origin of the atom-bearing region along the unit cell A, B, and C axes, in units of
  ///         spatial decomposition cells
  ///  4- 6   Extent of the atom-bearing region along the unit cell A, B, and C axes (cell indices
  ///         which exceed the system bounds will wrap), in units of spatial decomposition cells
  ///  7- 9   Origin of the grid-mapping region along the unit cell A, B, and C axes, in units of
  ///         spatial decomposition cells
  /// 10-12   Extent of the grid-mapping region along the unit cell A, B, and C axes, in units of
  ///         spatial decomposition cells
  /// 13-15   The number of spatial decomposition cells in the system along its unit cell A, B,
  ///         and C axes (copied from the associated CellGrid for convenience)
  ///  16     The total number of chains, the area of the cross section (pre-computed for
  ///         convenience)
  ///  17     The total size of the grid mapping region for this work unit (pre-computed for
  ///         convenience), a number of grid points (not spatial decomposition cells)
  /// 18-21   The lengths of the system's particle-mesh interaction grid along the A, B, and C unit
  ///         cell axes (indices 18-20), plus the offset for the system's grid in the concatenated
  ///         array.  This information is taken from the PMIGrid object's grid_dimensions array,
  ///         copied over to pack all such information into a single cache line.
  /// 22-24   Lengths of the system's spatial decomposition cell grid along the A, B, and C axes
  ///  25     The starting index of spatial decomposition cells in the system to which the work
  ///         unit pertains
  ///  26     The starting index of chains in the system to which the work unit pertains
  ///  27     The total number of warps assigned to each chain.
  ///  28     The total number of warp tasks.  This is the number of warps that will pile onto
  ///         each chain, times the total cross sectional area of the work unit.
  ///  29     The atom stride that each warp will take within a chain.
  Hybrid<uint> work_units;

  /// Pointer to the associated CellGrid object.  A PMIGrid cannot exist without a spatial
  /// decomposition, as the particle-mesh interaction grid is sized by the subdivision rate
  /// presented with the cell grid.  This pointer is recast from the original object.  The actual
  /// template types are stored in the member variables cg_tmat, cg_tacc, cg_tcalc, and cg_tcrd to
  /// recover the original pointer when needed.
  CellGrid<double, double, double, double4> *cg_pointer;

  // The following integers provide type ID numbers for the templated characteristics of the
  // attached CellGrid object.
  size_t cg_tmat;   ///< Data type for the local cell grid dimensions matrix (also a transformation
                    ///<   matrix taking fractional coordinates into Cartesian space, if needed in
                    ///<   such a way)
  size_t cg_tacc;   ///< Data type for accumulation of forces in the cell grid
  size_t cg_tcalc;  ///< Data type for transformation matrices taking Cartesian coordinates into
                    ///<   the coordinate space of one of the spatial decomposition cells
  size_t cg_tcrd;   ///< Data type for coordinates in the cell grid
  
  /// \brief Pointer to the associated coordinate synthesis
  PhaseSpaceSynthesis *poly_ps_pointer;

  /// \brief Validate either of the fixed-precision bit settings to ensure that it stays within the
  ///        allowed range.
  ///
  /// \param fp_bits  Number of bits after the decimal in the fixed-precision model
  void validateFixedPrecisionBits(int fp_bits) const;
  
  /// \brief Validate a system index within the associated synthesis.
  ///
  /// \param system_index  The system of interest
  void validateSystemIndex(int system_index) const;

  /// \brief Produce the maximum amount of memory that the __shared__ accumulation can use, based
  ///        on specifications of the available GPU.
  ///
  /// \param gpu  Details of the GPU that will run the calculations
  int findSharedBufferSize(const GpuDetails &gpu) const;
  
  /// \brief Compute the mapping halo in terms of spatial decomposition cells, for the current
  ///        order of B-spline interpolation and the attached cell grid.
  int computeMappingHalo() const;
  
  /// \brief Add a new work unit, as defined above, to the growing list.  The starting points of
  ///        the work unit's grid-mapping region, situated in the attached cell grid, are provided.
  ///        The complete extent of each work unit's atom-bearing region is calculated internally.
  ///
  /// \param result  The growing array of work units (modified and returned)
  /// \param sysid   Index of the system, within the attached synthesis, to which the work unit
  ///                applies
  /// \param pos_a   Position in the spatial decomposition cell grid at which the work unit begins
  ///                along the unit cell A axis
  /// \param pos_b   Position at which the work unit begins along the unit cell B axis
  /// \param pos_c   Position at which the work unit begins along the unit cell C axis
  /// \param gmap_a  Length of the work unit's grid-mapping region, measured in spatial
  ///                decomposition cells, along the unit cell's A axis
  /// \param gmap_b  Length of the work unit's grid-mapping region along the unit cell's B axis
  /// \param gmap_c  Length of the work unit's grid-mapping region along the unit cell's C axis
  /// \param cg_vc   A read-only abstract of the attached cell grid, needed for its non-templated
  ///                pointer to the spatial decomposition cell dimensions
  /// \param halo    The extent of the atom-mapping halo, the maximum number of spatial
  ///                decomposition cells which one might have to look down, back, or to the left
  ///                in order to find all atoms which might map their density to a grid point in
  ///                the current cell
  void addWorkUnit(std::vector<uint> *result, int sysid, int pos_a, int pos_b, int pos_c,
                   int gmap_a, int gmap_b, int gmap_c,
                   const CellGridReader<void, void, void, void> &cgr_vc, int halo = 1);

  /// \brief Check whether the topologies describing systems mapped in this object, the
  ///        fixed-precision bit count, and the non-bonded potential at hand will permit the use
  ///        of an abdriged accumulation approach which ignores overflow into the secondary 32-bit
  ///        accumulator.  To forego this check on the accumulation allows atomics to become "fire
  ///        and forget", as well as avoiding some arithmetic, which can enhance performance by
  ///        as much as 75%.
  void checkShortFormatViability();

  /// \brief Compute the number of grid points present in the largest work unit.
  void computeLargestWorkUnitGridPoints();
  
  /// \brief Unroll the inner data types of the attached cell grid and return the read-only
  ///        abstract for the attached CellGrid object.
  ///
  /// \param tier  Specify whether to obtain pointers, void-casted or otherwise, on the CPU host or
  ///              GPU device
  template <typename T, typename T4> const CellGridReader<void, void, void, void>
  unrollTemplateFreeCGReader(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
};

#ifdef STORMM_USE_HPC
/// \brief Launch the kernel to initialize the particle-mesh interaction grids (the object must be
///        configured for accumulation).
///
/// Overloaded:
///   - Provide the full GPU specifications
///   - Provide just a block count
///
/// \param pm_acc       Writeable, fixed-precision abstract of the PMIGrid object
/// \param gpu          Details of the available GPU, if initialization focuses on device memory
/// \param block_count  The number of blocks with which to launch initialization (useful if the
///                     full GPU specifications are not available)
/// \{
void launchPMIGridInitialization(PMIGridAccumulator *pm_acc, const GpuDetails &gpu);

void launchPMIGridInitialization(PMIGridAccumulator *pm_acc, int block_count);
/// \}

/// \brief Launch the kernel to convert the particle-mesh interaction grids from a fixed-precision
///        representation used in accumulation to a real-valued representation suitable for FFTs.
///        Overloading and descriptions of parameters follow from launchPMIGridInitialization()
///        above.
/// \{
void launchPMIGridRealConversion(PMIGridWriter *pm_wrt, const PMIGridAccumulator &pm_acc,
                                 const GpuDetails &gpu);

void launchPMIGridRealConversion(PMIGridWriter *pm_wrt, const PMIGridAccumulator &pm_acc,
                                 int block_count);
/// \}
#endif

} // namespace energy
} // namespace stormm

#include "pmigrid.tpp"

#endif
