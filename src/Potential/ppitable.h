// -*-c++-*-
#ifndef STORMM_PPI_TABLE_H
#define STORMM_PPI_TABLE_H

#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "Math/log_scale_spline.h"
#include "Math/math_enumerators.h"
#include "Parsing/parsing_enumerators.h"
#include "energy_enumerators.h"
#include "pme_util.h"

namespace stormm {
namespace energy {
  
using constants::ExceptionResponse;
using constants::PrecisionModel;
using symbols::amber_ancient_bioq;
using card::Hybrid;
using card::HybridKind;
using card::HybridTargetLevel;
using parse::NumberFormat;
using stmath::BasisFunctions;
using stmath::LogScaleSpline;
using stmath::LogSplineForm;
using stmath::LogSplineTable;
using stmath::TableIndexing;
using stmath::doublePrecisionSplineDetailMask;
using stmath::singlePrecisionSplineDetailMask;

/// \brief Abstract for the particle-particle interaction table, with pointers to energies and
///        forces and sizing constants for navigating the excluded versus non-excluded
///        interactions.  Like the LogScaleSpline abstract, this is read-only.
template <typename T, typename T4> struct PPIKit {

  /// \brief The constructor takes the customary list of pointers and critical constants.
  PPIKit(NonbondedTheme theme_in, BasisFunctions basis_in, TableIndexing lookup_in,
         int index_bound_in, int excl_offset_in, int index_shift_bits_in, ullint dp_detail_mask_in,
         uint sp_detail_mask_in, T arg_offset_in, const T4* energy_in, const T4* force_in,
         const T4* energy_excl_in, const T4* force_excl_in);

  /// \brief With all const members, the copy and move assignment operators are implicitly
  ///        deleted, but with no pointers to repair the copy and move constructors are legal.
  ///
  /// \param original  The original PPIKit tocopy or move
  /// \{
  PPIKit(const PPIKit &original) = default;
  PPIKit(PPIKit &&original) = default;
  /// \}

  const NonbondedTheme theme;   ///< The type of non-bonded potential and derivatives approximated
                                ///<   by the splines
  const BasisFunctions basis;   ///< The type of basis functions scaled by spline coefficients to
                                ///<   obtain the result
  const TableIndexing lookup;   ///< The method of transforming the argument to the underlying
                                ///<   function
  const int index_bound;        ///< The upper bound of spline indices in any given table.  This
                                ///<   is half the value of excl_offset.
  const int excl_offset;        ///< The offset to add to any spline index in the energy or force
                                ///<   arrays in order to obtain the corresponding excluded energy
                                ///<   or force spline coefficients
  const int index_shift_bits;   ///< The number of bits by which to shift a floating point number,
                                ///<   after interpreting it as an unsigned integer, to the right
                                ///<   before interpeting the result as a table index
  const ullint dp_detail_mask;  ///< A pre-computed bitmask useful for converting an argument of
                                ///<   the underlying benhmark function into a quantity for
                                ///<   double-precision splines to operate on
  const uint sp_detail_mask;    ///< A pre-computed bitmask useful for converting an argument of
                                ///<   the underlying benchmark function into a quantity for
                                ///<   single-precision splines to operate on
  const T arg_offset;           ///< Offset used in transforming arguments of the underlying
                                ///<   benchmark function into spline table indices
  const T4* energy;             ///< Table of energy coefficients.  Add excl_offset to any index
                                ///<   for the raw energy to get the corresponding index in the
                                ///<   excluded energy table.  Add index_bound to any index for the
                                ///<   raw energy to get the raw force.
  const T4* force;              ///< Table of force coefficients.  Add excl_offset to any index
                                ///<   for the raw force to get the corresponding index in the
                                ///<   excluded force table.
  const T4* energy_excl;        ///< Table of excluded energy coefficients
  const T4* force_excl;         ///< Table of excluded force coefficients
};

/// \brief An equivalent abstract which delivers the tabulated data with all tuples broken into
///        separate arrays of their individual components.  This abstract is obtained with distinct
///        accessors from PPITable, below.  The "e" is for "element-wise."
template <typename T> struct PPIeKit {

  /// \brief The constructor takes the customary list of pointers and critical constants.
  PPIeKit(NonbondedTheme theme_in, BasisFunctions basis_in, TableIndexing lookup_in,
          int index_bound_in, int excl_offset_in, int index_shift_bits_in,
          ullint dp_detail_mask_in, uint sp_detail_mask_in, T arg_offset_in, const T* energy_x_in,
          const T* energy_y_in, const T* energy_z_in, const T* energy_w_in, const T* force_x_in,
          const T* force_y_in, const T* force_z_in, const T* force_w_in, const T* energy_excl_x_in,
          const T* energy_excl_y_in, const T* energy_excl_z_in, const T* energy_excl_w_in,
          const T* force_excl_x_in, const T* force_excl_y_in, const T* force_excl_z_in,
          const T* force_excl_w_in);

  /// \brief With all const members, the copy and move assignment operators are implicitly
  ///        deleted, but with no pointers to repair the copy and move constructors are legal.
  ///
  /// \param original  The original PPIKit tocopy or move
  /// \{
  PPIeKit(const PPIeKit &original) = default;
  PPIeKit(PPIeKit &&original) = default;
  /// \}

  const NonbondedTheme theme;   ///< The type of non-bonded potential and derivatives approximated
                                ///<   by the splines
  const BasisFunctions basis;   ///< The type of basis functions scaled by spline coefficients to
                                ///<   obtain the result
  const TableIndexing lookup;   ///< The method of transforming the argument to the underlying
                                ///<   function
  const int index_bound;        ///< The upper bound of spline indices in any given table.  This
                                ///<   is half the value of excl_offset.
  const int excl_offset;        ///< The offset to add to any spline index in the energy or force
                                ///<   arrays in order to obtain the corresponding excluded energy
                                ///<   or force spline coefficients
  const int index_shift_bits;   ///< The number of bits by which to shift a floating point number,
                                ///<   after interpreting it as an unsigned integer, to the right
                                ///<   before interpeting the result as a table index
  const ullint dp_detail_mask;  ///< A pre-computed bitmask useful for converting an argument of
                                ///<   the underlying benhmark function into a quantity for
                                ///<   double-precision splines to operate on
  const uint sp_detail_mask;    ///< A pre-computed bitmask useful for converting an argument of
                                ///<   the underlying benchmark function into a quantity for
                                ///<   single-precision splines to operate on
  const T arg_offset;           ///< Offset used in transforming arguments of the underlying
                                ///<   benchmark function into spline table indices
  const T* energy_x;            ///< Table of energy function spline x coefficients.  Add
                                ///<   excl_offset to any index for the raw energy to get the
                                ///<   corresponding index for the x coefficient in the excluded
                                ///<   energy table.
  const T* energy_y;            ///< Table of energy function spline y coefficients
  const T* energy_z;            ///< Table of energy function spline z coefficients
  const T* energy_w;            ///< Table of energy function spline w coefficients
  const T* force_x;             ///< Table of force function spline x coefficients.  Add
                                ///<   excl_offset to obtain the force x coefficient for an
                                ///<   excluded interaction.
  const T* force_y;             ///< Table of force function spline y coefficients
  const T* force_z;             ///< Table of force function spline z coefficients
  const T* force_w;             ///< Table of force function spline w coefficients
  const T* energy_excl_x;       ///< Excluded interation energy function spline x coefficients
  const T* energy_excl_y;       ///< Excluded interation energy function spline y coefficients
  const T* energy_excl_z;       ///< Excluded interation energy function spline z coefficients
  const T* energy_excl_w;       ///< Excluded interation energy function spline w coefficients
  const T* force_excl_x;        ///< Excluded interation energy function spline x coefficients
  const T* force_excl_y;        ///< Excluded interation energy function spline y coefficients
  const T* force_excl_z;        ///< Excluded interation energy function spline z coefficients
  const T* force_excl_w;        ///< Excluded interation energy function spline w coefficients
};
  
/// \brief A tabulated non-bonded potential, with or without exclusions, to be used in the context
///        of particle-particle, particle-mesh calculations.  The key is to create two tables, one
///        for the non-excluded form of the interaction and the other for the excluded form.  The
///        entries for each table will then be concatenated, such that all non-excluded
///        interactions are contiguous and then all exclude interactions are contiguous.  The
///        offset for accessing an excluded interaction based on an index calculated from a
///        particle-particle distance is stored alongside the tabulated splines in the abstract.
class PPITable {
public:

  /// \brief The constuctor can accept all of the arguments that might be useful for making a
  ///        LosScaleSpline, or a LogScaleSpline itself.  Tables for both the otential and the
  ///        derivative will be computed.
  /// \{
  PPITable(NonbondedTheme theme_in = NonbondedTheme::ELECTROSTATIC,
           BasisFunctions basis_set_in = BasisFunctions::POLYNOMIAL,
           TableIndexing indexing_method_in = TableIndexing::SQUARED_ARG,
           double cutoff_in = default_pme_cutoff, double argument_offset_in = 0.0,
           double dsum_tol_in = default_dsum_tol, int mantissa_bits_in = 5,
           double coulomb_in = amber_ancient_bioq, double min_range_in = 0.015625);
  
  template <typename T4>
  PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b,
           const LogScaleSpline<T4> &spl_c, const LogScaleSpline<T4> &spl_d,
           double cutoff_in = default_pme_cutoff);

  template <typename T4>
  PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b,
           const LogScaleSpline<T4> &spl_c);

  template <typename T4>
  PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b);

  template <typename T4>
  PPITable(const LogScaleSpline<T4> &spl_a);
  /// \}
  
  /// \brief The presence of POINTER-kind Hybrid objects implies pointer repairs that require the
  ///        copy and move constructors as well as assignment operators to be spelled out.
  ///
  /// \param original  THe original object to copy or move
  /// \param other     Another object placed on the right hand side of the assignment statement
  /// \{
  PPITable(const PPITable &original);
  PPITable(PPITable &&original);
  PPITable& operator=(const PPITable &original);
  PPITable& operator=(PPITable &&original);
  /// \}

  /// \brief Get the non-bonded potential described by an object of this class.
  NonbondedTheme getTheme() const;

  /// \brief Get the cutoff on particle-particle interactions.
  double getCutoff() const;

  /// \brief Get the maximum range of the spline.
  double getMaximumRange() const;

  /// \brief Get the indexing argument offset, consistent across all tables.
  double getIndexingOffset() const;
  
  /// \brief Get the direct sum tolerance.
  double getDirectSumTolerance() const;

  /// \brief Get the number of bits of the mantissa used for table indexing.
  int getBitStride() const;

  /// \brief Get the Ewald coefficient used to perform the switching between short- and
  ///        long-ranged potentials.
  double getEwaldCoefficient() const;

  /// \brief Get the Gaussian RMS sigma parameter (one half the inverse of the Ewald coefficient)
  ///        used to perform the switching between short- and long-ranged potentials.
  double getGaussianWidth() const;

  /// \brief Get the spline index that would be accessed in order to evaluate a given argument.
  ///
  /// \param arg  The value of the spline table indexing argument
  int getTableIndex(double arg, PrecisionModel prec = PrecisionModel::SINGLE) const;

  /// \brief Get the spline index that would be accessed in order to evaluate a given argument.
  ///
  /// \param arg  The value of the spline table indexing argument
  int getTableIndexByRealArg(double arg, PrecisionModel prec = PrecisionModel::SINGLE) const;

  /// \brief Evaluate the spline table at a given inter-particle distance expressed in the manner
  ///        used for the class object's table indexing.  For example, if the object indexes its
  ///        tables by the squared value of the displacement, enter 49.0 A^2 in order to get the
  ///        force, energy, or excluded form thereof for two particles separated by a distance of
  ///        7.0 A.
  ///
  /// \param arg                   The value of the spline table indexing argument
  /// \param kind                  Indicate whether to evaluate the potential, the force, or an
  ///                              excluded form of either quantity
  /// \param prec                  The precision of the table from which to retrieve calculations
  ///                              (also implies a precision model for the spline calculation)
  /// \param use_elemental_tables  Indicate whether to use elementwise tables, in which each of the
  ///                              spline coefficients are accessed individually, or tuple tables
  ///                              in which they are all accessed at once.  The default behavior is
  ///                              to use tuple-baed tables.
  double evaluate(double arg, LogSplineForm kind, PrecisionModel prec = PrecisionModel::SINGLE,
                  bool use_elemental_tables = false) const;

  /// \brief Evaluate the spline table at a given inter-particle distance expressed as the absolute
  ///         value of the displacement between two particles.  Descriptions of input parameters
  ///         otherwise follow from evaluate(), above.
  double evaluateByRealArg(double arg, LogSplineForm kind,
                           PrecisionModel prec = PrecisionModel::SINGLE,
                           bool use_elemental_tables = false) const;
  
  /// \brief Get the double-precision abstract for use of the splines in a C programming style.
  ///
  /// \param tier  Indicate whether to obtain pointers for data on the CPU host or GPU device
  const PPIKit<double, double4> dpData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the single-precision abstract for use of the splines in a C programming style.
  ///
  /// \param tier  Indicate whether to obtain pointers for data on the CPU host or GPU device
  const PPIKit<float, float4> spData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the double-precision abstract for use of the splines in a C programming style,
  ///        with all spline tuples broken into their individual components.  This abstract may be
  ///        more suitable for kernels where register pressure constrains performance.
  ///
  /// \param tier  Indicate whether to obtain pointers for data on the CPU host or GPU device
  const PPIeKit<double> dpeData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the single-precision abstract for use of the splines in a C programming style,
  ///        with all spline tuples broken into their individual components.  This abstract may be
  ///        more suitable for kernels where register pressure constrains performance.
  ///
  /// \param tier  Indicate whether to obtain pointers for data on the CPU host or GPU device
  const PPIeKit<float> speData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

#ifdef STORMM_USE_HPC
  /// \brief Upload data to the GPU device.
  void upload();
  
  /// \brief Download data from the GPU device.
  void download();
#endif

private:

  NonbondedTheme theme;           ///< The type of non-bonded potential encoded in the tabulated
                                  ///<   coefficients.
  BasisFunctions basis_set;       ///< The type of basis functions to use in spline approximations
  TableIndexing indexing_method;  ///< The indexing method that converts the argument of the
                                  ///<   underlying function into a spline table element
  double cutoff;                  ///< Maximum inter-particle distance for which interactions are
                                  ///<   valid
  double max_range;               ///< The maximum interparticle distance for which the table can
                                  ///<   produce a result
  double min_range;               ///< The minimum range for spline computations.  Below this
                                  ///<   range, the spline values will be assumed to remain
                                  ///<   constant at the lowest explicitly computed value.
  double argument_offset;         ///< Offset that will be applied to the indexing argument when
                                  ///<   accessing spline tables, e.g. if the table is indexed by
                                  ///<   r2 plus 0.5, add 0.5 to the square of the inter-particle
                                  ///<   distance before extracting the table index
  double dsum_tol;                ///< The direct sum tolerance, a proportion of each interaction
                                  ///<   that is discarded beginning at the particle-particle pair
                                  ///<   cutoff
  double ew_coeff;                ///< The Ewald coefficient that governs the potential splitting
                                  ///<   between short- and long-ranged components
  int mantissa_bits;              ///< Number of mantissa bits used to index into various tables
  double coulomb;                 ///< The object's take on Coulomb's constant
  int dp_exclusion_offset;        ///< The offset between excluded and non-excluded spline table
                                  ///<   segments for DOUBLE-precision splines, when computing
                                  ///<   either energies or forces.
  int sp_exclusion_offset;        ///< The offset between excluded and non-excluded spline table
                                  ///<   segments for SINGLE-precision splines, when computing
                                  ///<   either energies or forces.
  ullint dp_detail_bitmask;       ///< Detail mask for extracting the non-indexing bits from a
                                  ///<   spline table based on DOUBLE precision float64_t numbers.
  uint sp_detail_bitmask;         ///< Detail mask for extracting the non-indexing bits from a
                                  ///<   spline table based on SINGLE precision float32_t numbers.

  // The main arrays, in full (double) precision
  Hybrid<double4> energy;  ///< Double-precision coefficients for the energy, without exclusions
  Hybrid<double4> force;   ///< Double-precision coefficients for the force, without exclusions

  /// \brief Double-precision coefficients for the energy, with excluded 1:4 interactions
  Hybrid<double4> energy_with_exclusions;

  /// \brief Double-precision coefficients for the force, with excluded 1:4 interactions
  Hybrid<double4> force_with_exclusions;

  // The main arrays, in single precision
  Hybrid<float4> sp_energy;  ///< Single-precision coefficients for the energy, without exclusions
  Hybrid<float4> sp_force;   ///< Single-precision coefficients for the force, without exclusions
  
  /// Single-precision coefficients for the energy, with excluded 1:4 interactions
  Hybrid<float4> sp_energy_with_exclusions;

  /// Single-precision coefficients for the force, with excluded 1:4 interactions
  Hybrid<float4> sp_force_with_exclusions;

  /// Coefficients of the tables computed in double precision.  This is an ARRAY-kind Hybrid
  /// targeted by the double-precision tuple POINTER-kind Hybrid arrays above.
  Hybrid<double4> coeffs;

  /// Coefficients of the tables computed in single precision.  This is an ARRAY-kind Hybrid
  /// targeted by the single-precision tuple POINTER-kind Hybrid arrays above.
  Hybrid<float4> sp_coeffs;

  /// Double-precision coefficients for the scalar, first, second, and third-order portions of
  /// each cubic spline function, without exclusions.  These are the "x", "y", "z", and "w"
  /// components of the energy member variable above.
  /// \{
  Hybrid<double> energy_x;
  Hybrid<double> energy_y;
  Hybrid<double> energy_z;
  Hybrid<double> energy_w;
  /// \}

  /// Double-precision coefficients for the scalar, first, second, and third-order portions of
  /// each cubic spline function, without exclusions.  These are the "x", "y", "z", and "w"
  /// components of the energy member variable above.
  /// \{
  Hybrid<double> force_x;
  Hybrid<double> force_y;
  Hybrid<double> force_z;
  Hybrid<double> force_w;
  /// \}

  /// Double-precision coefficients for the scalar, first, second, and third-order portions of
  /// each excluded energy cubic spline function, without exclusions.  These are the "x", "y", "z",
  /// and "w" components of the energy_with_exclusions member variable above.
  /// \{
  Hybrid<double> energy_with_excl_x;
  Hybrid<double> energy_with_excl_y;
  Hybrid<double> energy_with_excl_z;
  Hybrid<double> energy_with_excl_w;
  /// \}

  /// Double-precision coefficients for the scalar, first, second, and third-order portions of
  /// each excluded force cubic spline function, without exclusions.  These are the "x", "y", "z",
  /// and "w" components of the force_with_exclusions member variable above.
  /// \{
  Hybrid<double> force_with_excl_x;
  Hybrid<double> force_with_excl_y;
  Hybrid<double> force_with_excl_z;
  Hybrid<double> force_with_excl_w;
  /// \}

  /// Double-precision coefficients for the scalar, first, second, and third-order portions of
  /// each cubic spline function, without exclusions.  These are the "x", "y", "z", and "w"
  /// components of the energy member variable above.
  /// \{
  Hybrid<float> sp_energy_x;
  Hybrid<float> sp_energy_y;
  Hybrid<float> sp_energy_z;
  Hybrid<float> sp_energy_w;
  /// \}

  /// Double-precision coefficients for the scalar, first, second, and third-order portions of
  /// each cubic spline function, without exclusions.  These are the "x", "y", "z", and "w"
  /// components of the energy member variable above.
  /// \{
  Hybrid<float> sp_force_x;
  Hybrid<float> sp_force_y;
  Hybrid<float> sp_force_z;
  Hybrid<float> sp_force_w;
  /// \}

  /// Double-precision coefficients for the scalar, first, second, and third-order portions of
  /// each excluded energy cubic spline function, without exclusions.  These are the "x", "y", "z",
  /// and "w" components of the energy_with_exclusions member variable above.
  /// \{
  Hybrid<float> sp_energy_with_excl_x;
  Hybrid<float> sp_energy_with_excl_y;
  Hybrid<float> sp_energy_with_excl_z;
  Hybrid<float> sp_energy_with_excl_w;
  /// \}

  /// Double-precision coefficients for the scalar, first, second, and third-order portions of
  /// each excluded force cubic spline function, without exclusions.  These are the "x", "y", "z",
  /// and "w" components of the force_with_exclusions member variable above.
  /// \{
  Hybrid<float> sp_force_with_excl_x;
  Hybrid<float> sp_force_with_excl_y;
  Hybrid<float> sp_force_with_excl_z;
  Hybrid<float> sp_force_with_excl_w;
  /// \}

  /// Coefficients of all tables computed in double-precision, with all tuples broken into their
  /// component members.  This is an ARRAY-kind Hybrid targeted by the double-precision tuple
  /// POINTER-kind Hybrid arrays above.
  Hybrid<double> elemental_coeffs;
  
  /// Coefficients of all tables computed in single-precision, with all tuples broken into their
  /// component members.  This is, again, an ARRAY-kind Hybrid targeted by the double-precision
  /// tuple POINTER-kind Hybrid arrays above.
  Hybrid<float> sp_elemental_coeffs;
  
  /// \brief Determine the theme of the particle-particle interaction table based on the form of
  ///        one of the input spline tables.
  ///
  /// \param spl  The input spline table, assumed to be representative of the forms of other spline
  ///             tables fed to the PPITable object
  template <typename T4> NonbondedTheme findTheme(const LogScaleSpline<T4> &spl) const;

  /// \brief Build the tables for energy and force computations from pairwise interactions and pack
  ///        the results into a vector for later use.
  template <typename T4> std::vector<LogScaleSpline<T4>> buildAllSplineTables() const;
  
  /// \brief Find a potential function, without exclusions, among the spline tables provided.
  ///
  /// \param spl_a  The first spline table provided
  /// \param spl_b  The second spline table provided
  /// \param spl_c  The third spline table provided
  /// \param spl_d  The fourth spline table provided
  template <typename T4>
  const LogScaleSpline<T4>& findNonExclPotential(const LogScaleSpline<T4> &spl_a,
                                                 const LogScaleSpline<T4> &spl_b,
                                                 const LogScaleSpline<T4> &spl_c,
                                                 const LogScaleSpline<T4> &spl_d) const;

  /// \brief Find a potential function, with non-bonded 1:4 exclusions, among the spline tables
  ///        provided.  Descriptions of input parameters follow from findNonExclPotential() above.
  template <typename T4>
  const LogScaleSpline<T4>& findExclPotential(const LogScaleSpline<T4> &spl_a,
                                              const LogScaleSpline<T4> &spl_b,
                                              const LogScaleSpline<T4> &spl_c,
                                              const LogScaleSpline<T4> &spl_d) const;

  /// \brief Find a force function, without exclusions, among the spline tables provided.
  ///        Descriptions of input parameters follow from findNonExclPotential() above.
  template <typename T4>
  const LogScaleSpline<T4>& findNonExclForce(const LogScaleSpline<T4> &spl_a,
                                             const LogScaleSpline<T4> &spl_b,
                                             const LogScaleSpline<T4> &spl_c,
                                             const LogScaleSpline<T4> &spl_d) const;

  /// \brief Find a force function, with non-bonded 1:4 exclusions, among the spline tables
  ///        provided.  Descriptions of input parameters follow from findNonExclPotential() above.
  template <typename T4>
  const LogScaleSpline<T4>& findExclForce(const LogScaleSpline<T4> &spl_a,
                                          const LogScaleSpline<T4> &spl_b,
                                          const LogScaleSpline<T4> &spl_c,
                                          const LogScaleSpline<T4> &spl_d) const;

  /// \brief Check a bit in a mask to signify the presence of a particular, necessary function.
  ///
  /// \param spl_x  A spline containing one function relevant to the PPITable
  template <typename T4> uint checkPriority(const LogScaleSpline<T4> &spl_x) const;

  /// \brief Check whether two splines share the same length and density parameters such that they
  ///        would be compatible to place in a table together.
  ///
  /// \param spl_a  The first spline
  /// \param spl_b  The second spline to compare against the first
  template <typename T4> void checkSplineCompatibility(const LogScaleSpline<T4> &spl_a,
                                                       const LogScaleSpline<T4> &spl_b) const;
  
  /// \brief Determine the highest-priority missing functional form.
  ///
  /// \param holdings  A bitmask checked 1 for the availability of relevant functions: 0x1 for
  ///                  the potential, 0x2 for the excluded potential, 0x4 for the force, 0x8 for
  ///                  the excluded force
  LogSplineForm findMissingForm(const uint holdings) const;
  
  /// \brief Construct the necessary splines for this aggregate table, based on one to three
  ///        other splines which specify the type.  The priority system is first to find the
  ///        non-excluded energy, then the excluded energy, then the non-excluded force, and
  ///        finally the excluded force (between particles).  Various overloads work given up to
  ///        three available logarithmic spline tables.
  ///
  /// \param spl_a  The first of up to three known spline tables
  /// \param spl_b  The second of up to three known spline tables
  /// \param spl_c  The third and last known spline table
  /// \{
  template <typename T4>
  LogScaleSpline<T4> getTablePriority(const LogScaleSpline<T4> &spl_a,
                                      const LogScaleSpline<T4> &spl_b,
                                      const LogScaleSpline<T4> &spl_c) const;

  template <typename T4>
  LogScaleSpline<T4> getTablePriority(const LogScaleSpline<T4> &spl_a,
                                      const LogScaleSpline<T4> &spl_b) const;

  template <typename T4>
  LogScaleSpline<T4> getTablePriority(const LogScaleSpline<T4> &spl_a) const;
  /// \}

  /// \brief Fill the coefficients table with the energy, force, excluded energy, and excluded
  ///        force splines.
  ///
  /// \param u    The energy spline table
  /// \param du   The force spline table
  /// \param ux   The excluded energy spline table
  /// \param dux  The excluded force spline table
  template <typename T4>
  void populateCoefficients(const LogScaleSpline<T4> &u, const LogScaleSpline<T4> &du,
                            const LogScaleSpline<T4> &ux, const LogScaleSpline<T4> &dux);
};

} // namespace energy
} // namespace stormm

#include "ppitable.tpp"

#endif
