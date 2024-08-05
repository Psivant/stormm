// -*-c++-*-
#ifndef STORMM_LOG_SCALE_SPLINE_H
#define STORMM_LOG_SCALE_SPLINE_H

#include <cmath>
#include <vector>
#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/matrix_ops.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "Potential/pme_util.h"
#include "Reporting/error_format.h"
#include "math_enumerators.h"

namespace stormm {
namespace stmath {

using card::GpuDetails;
using card::Hybrid;
using card::HybridTargetLevel;
using constants::ExceptionResponse;
using constants::PrecisionModel;
using symbols::amber_ancient_bioq;
using energy::minimum_ewald_coefficient;
using energy::maximum_ewald_coefficient;
using parse::minimalRealFormat;
using parse::NumberFormat;
using parse::realToString;

/// \brief Default settings for the logarthmic spline tables.  The range limits are given in units
///        of Angstroms, not Angstroms^2.
/// \{
constexpr int default_logtab_mantissa_bits = 5;
constexpr double default_logtab_max_range = 64.0;
constexpr double default_logtab_min_range = 0.00390625;
/// \}
  
/// \brief Abstract for the logarithmic splined function object, containing the table and critical
///        constants for evaluating it.
template <typename T4> struct LogSplineTable {

  /// \brief As with other abstracts, the constructor takes a straight list of inputs for each
  ///        member variable.
  LogSplineTable(BasisFunctions basis_in, TableIndexing lookup_in, int detail_bits_in,
                 int index_bound_in, uint sp_detail_mask_in, ullint dp_detail_mask_in,
                 const float arg_offset_in, const T4* table_in);

  /// \brief The const-ness of member variables implicitly deletes the copy and move assignment
  ///        operators, but the default copy and move constructors are valid.
  /// \{
  LogSplineTable(const LogSplineTable &original) = default;
  LogSplineTable(LogSplineTable &&original) = default;
  /// \}

  /// The type of basis functions f(x), g(x), h(x), and v(x) in the expression evaluated by the
  /// spline:
  ///
  /// U(x) = A f(x) + B g(x) + C h(x) + D v(x)
  const BasisFunctions basis;

  /// The method of transforming the argument of the underlying benchmark function into a table
  /// index for selecting the correct spline
  const TableIndexing lookup;
  
  /// The number of bits of the corresponding floating point format of the range argument.  For a
  /// double4 LogSplineTable, the range argument is double and the number of detail bits (on any
  /// architecture STORMM is prepared for use on) should be 52 minus the number of mantissa bits
  /// used in the index.  For float4 data and float range arguments, the number of index bits will
  /// be 23 minus the number of mantissa bits.
  const int detail_bits;

  /// The maximum index of the table, for bounds checking
  const int index_bound;

  /// The following detail masks guide the conversion of a range argument for the function of
  /// interest into an argument for polynomial basis functions.  Both single- and double-precision
  /// forms of these masks must be included, as they cannot be cast to void (only pointers may be
  /// stripped of their templated characteristics).  The detail masks are bitmasks with the lowest
  /// detail_bits set to one and all other bits set to zero.  In effect, the range argument becomes
  /// divided into two parts, the high bits being an index and the low bits being the detail left
  /// over to submit to the interpolant's basis functions.
  ///
  /// When the underlying LogScaleSpline object works on the BasisFunctions::POLYNOMIAL setting,
  /// the process is to take the range argument for the function of interest (e.g. the squared
  /// distance between two particles for computing the derivative of the Ewald direct space sum),
  /// mask off the final bits into an unsigned integer of equal size using the appropriate detail
  /// mask, shift these bits forward by the number of mantissa bits used to determine the table
  /// index (to put the remaining detail at the front of the mantissa, protecting it from roundoff
  /// error when computing the square or cube term), and replace the exponent bits with a mask
  /// indicating the range [1, 2).  For IEEE float, the exponent bits are 0x3f800000 and for IEEE
  /// double the exponent bits are 0x3ff000000000000ULL (placing 127 and 1023 in the exponents of
  /// each number, respectively).
  /// \{
  const uint sp_detail_mask;
  const ullint dp_detail_mask;
  /// \}

  /// The offset added to whatever argument (squared or unsquared value of the argument to the
  /// benchmark function) indexes into the table.  Applied only if the indexing method is
  /// ARG_OFFSET or SQ_ARG_OFFSET.  Represented as a float, as it must be a small power of two and
  /// thus entails no loss of information to represent it in the shorter format.
  const float arg_offset;
  
  /// The table of cubic spline coefficients
  const T4* table;
};
  
/// \brief A logarithmic spline can be very useful for interpolating functions that are steepest
///        at low values of the argument, provided that the argument is never negative.  This is
///        a useful class of functions for describing inter-atomic potentials.  The strategy is
///        to utilize the logarithmic nature (and IEEE format) of a floating point number, taking
///        the highest N bits of the mantissa and the exponent as an unsigned integer (because
///        measurements of absolute distance will never be negative, the highest bit (the sign bit)
///        will always be set to zero).  It is most useful if the square root operation can be
///        avoided altogether, indexing the table by the square of the function argument.  This
///        creates a more accurate table as well as saving arithmetic, although obtaining function
///        derivatives (as a function of the argument, not its square) from the table becomes more
///        expensive.  In fact, this object abandons continuous derivatives in the tables and
///        instead splines the derivative of the desired function as a separate logarithmic table.
///
///        The process generates 2^N spline segments in the ranges [0, 2^-k), [2^-k, 2^(1-k)),
///        [2^(1-k), 2^(2-k)), and so on for all k representable by floating point format's
///        exponent.  For N=5, there are 32 spline segments covering the ranges [1/4, 1/2),
///        [1/2, 1), [1, 2), [2, 4), and so on.  While a great deal of the table is spent covering
///        infinitesimal increments over very negative powers of 2, the quantity of such numbers is
///        limited by the size of the exponent (256 for 32-bit floats and 2048 for 64-bit floats in
///        IEEE format), and these values are hardly ever accessed.  For N = 5 in 32-bit floats,
///        the cubic spline coefficients occupy 128kB of space, of which splines to cover the range
///        [0.5, 16.0) occupy 2.5kB.  For N = 8 in 64-bit floats, the cubic spline coefficients
///        occupy 8MB of space and splines to cover the range [0.5, 16.0) occupy 40kB.  The latter
///        may be of use on CPU resources with higher cache, or in GPU contexts where caching can
///        be well managed, but the former is quite useful in many GPU cases.
template <typename T4> class LogScaleSpline {
public:

  /// \brief The constructor needs an indication of the functional form to render.  The number of
  ///        bits will be decided by the data type if not indicated by the developer.
  ///        Interpolation is fixed at cubic, for simplicity and speed.  Additional input may be
  ///        used to specify a particular function to spline.
  ///
  /// \param custom_form_in  If provided, this collection of arrays will provide the known values
  ///                        of the function to be splined (in each tuple's "y" member) and the
  ///                        value of the function argument for which the function is known (in
  ///                        each tuple's "x" member)
  /// \param max_range_in    This will translate into the maximum_range member variable, after some
  ///                        checking and refinement
  /// \param min_range_in    This will translate into the minimum_absolute_range member variable,
  ///                        after some checking and refinement
  /// \param policy          Indicate what to do in the event of bad input, i.e. if custom spline
  ///                        data is not able to fulfill the range requested
  /// \{
  LogScaleSpline(LogSplineForm target_form_in, double ewald_coefficient_in,
                 double coulomb_constant_in = amber_ancient_bioq,
                 int mantissa_bits_in = default_logtab_mantissa_bits,
                 double max_range_in = default_logtab_max_range,
                 double min_range_in = default_logtab_min_range,
                 TableIndexing indexing_method_in = TableIndexing::SQUARED_ARG,
                 BasisFunctions basis_set_in = BasisFunctions::POLYNOMIAL,
                 int ulp_optimization_depth_in = 2,
                 float indexing_offset_in = 0.0,
                 ExceptionResponse policy = ExceptionResponse::DIE);

  LogScaleSpline(LogSplineForm target_form_in,
                 const std::vector<std::vector<double2>> &custom_form_in,
                 int mantissa_bits_in = default_logtab_mantissa_bits,
                 double max_range_in = default_logtab_max_range,
                 double min_range_in = default_logtab_min_range,
                 TableIndexing indexing_method_in = TableIndexing::SQUARED_ARG,
                 BasisFunctions basis_set_in = BasisFunctions::POLYNOMIAL,
                 int ulp_optimization_depth_in = 2,
                 float indexing_offset_in = 0.0,
                 ExceptionResponse policy = ExceptionResponse::WARN);
  /// \}

  /// \brief With no pointers to repair and no const members, the default copy and move
  ///        constructors as well as assignment operators are all valid.
  ///
  /// \param original  The existing object to copy or move
  /// \param other     Another object placed on the right hand side of the assignment statement
  /// \{
  LogScaleSpline(const LogScaleSpline &original) = default;
  LogScaleSpline(LogScaleSpline &&original) = default;
  LogScaleSpline& operator=(const LogScaleSpline &original) = default;
  LogScaleSpline& operator=(LogScaleSpline &&original) = default;
  /// \}

  /// \brief Get the target form of the function represented by the table.
  LogSplineForm getForm() const;

  /// \brief Get the method of indexing the table, an expression based on the target function's
  ///        argument, e.g. the squared distance.
  TableIndexing getIndexingMethod() const;

  /// \brief Get the basis functions used to construt each spline.
  BasisFunctions getBasisSet() const;
  
  /// \brief Get the number of leading bits in the mantissa used to determine the array index.
  int getBitStride() const;

  /// \brief Get the number of spline segments per stride.  This is two raised to the power of the
  ///        result of getBitStride().
  int getSplineDensity() const;

  /// \brief Get the Ewald coefficient in use for this spline table.  A check will be made to
  ///        ensure that the spline is of the proper type to have a valid Ewald coefficient.
  double getEwaldCoefficient() const;

  /// \brief Get the value of Coulomb's constant in use by the spline.  A check will be made to
  ///        ensure that the spline is of the proper type to have a valid Coulomb constant.
  double getCoulombConstant() const;
  
  /// \brief Get the maximum range of the spline.
  double getMaximumRange() const;

  /// \brief Get the minimum absolute range of the spline.
  double getMinimumRange() const;

  /// \brief Get the indexing offset.
  float getIndexingOffset() const;

  /// \brief Get the optimization depth.
  int getOptimizationDepth() const;
  
  /// \brief Get the spline index for a particular value of the argument of the tabulated function
  ///        which has already been transformed to meet the object's indexing system.
  ///
  /// \param r  The argument to the tabulated function.  The exact meaning will be interpreted
  ///           based on the table's internal indexing_method setting.
  int getSplineIndex(double r) const;

  /// \brief Get the spline index for a particular value of the argument of the tabulated function.
  ///        Transformation of the input will occur internally based on how the object indexes its
  ///        internal spline table, in a manner consistent with the object's precision model.
  ///
  /// \param r  The argument to the tabulated function, expressed on its natural axis
  int getSplineIndexByRealArg(double r) const;

  /// \brief Get the estimated error as a function of the distance at which the splines are
  ///        evaluated.  Mean error is returned in the "x" member of the tuple while standard
  ///        deviation about the mean error is given in the "y" member of the tuple.  The maximum
  ///        and error value is returned in the tuple's "z" member.
  ///
  /// Overloaded:
  ///   - Get the overall error estimate
  ///   - Get the error estimate for the spline segment comprising a particular value of distance
  ///
  /// \param r  The distance of interest
  /// \{
  double3 getErrorEstimate() const;
  double3 getErrorEstimate(double r) const;
  /// \}

  /// \brief Get the value of the underlying function based on a value of the function argument.
  ///        This accessor will throw a warning or an error if trying to access the a table built
  ///        to be accessed using the squared value of the function.
  ///
  /// \param r  The value of at which to evaluate the underlying function
  double interpolateByValue(double r) const;

  /// \brief Get the value of the underlying function based on a squared value of the function
  ///        argument.  This accessor will throw a warning or an error if trying to access the a
  ///        table built to be accessed using the value of the function.
  ///
  /// \param r2  The squared value at which to evaluate the underlying function
  double interpolateBySquaredValue(double r2) const;

  /// \brief Evaluate the underlying function represented by the spline at a point expressed in
  ///        the object's native indexing system (meaning, if the object indexes by the squared
  ///        argument of the underlying function, then the input will reflect r^2 rather than r).
  ///
  /// \param r  Argument to the underlying function, transformed as appropriate for indexing
  double evaluate(double r) const;

  /// \brief Evaluate the underlying function represented by the spline at a point on the natural
  ///        axis of the underlying function.  Transformation of this argument will be handled
  ///        internally, and in a manner that reflects the object's precision model, before
  ///        returning the result.
  ///
  /// \param r  Argument to the underlying function
  double evaluateByRealArg(double r) const;  
  
  /// \brief Get the object's abstract for the CPU host or GPU device in either of two precision
  ///        modes.
  ///
  /// \brief tier  Indicate whether to retrieve pointers on the GPU device or CPU host
  const LogSplineTable<T4> data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the object's abstract with templating removed.
  const LogSplineTable<void>
  templateFreeData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a pointer to the original object.
  const LogScaleSpline<T4>* getSelfPointer() const;
  
#ifdef STORMM_USE_HPC
  /// \brief Upload the table data to the device.
  void upload();
  
  /// \brief Download the table data from the device.
  void download();
#endif
  /// \brief Set the minimum significant range.  This is purely for error analysis and will not
  ///        affect the array indexing or other details of the spline.  The value is set to 0.5
  ///        at initialization.
  ///
  /// \param min_range_in  The minimum range to set
  void setMinimumSignificantRange(double min_range_in);

private:
  LogSplineForm target_form;                ///< The target functional form, which can be a number
                                            ///<   of known forms of a custom potential
  TableIndexing indexing_method;            ///< The manner in which the table will be indexed,
                                            ///<   indicating some power of the function argument
  BasisFunctions basis_set;                 ///< Indicate the form of the basis functions for the
                                            ///<   interpolant.  Any code utilizing the spline will
                                            ///<   need to implement these basis functions, and
                                            ///<   one set may be more convenient than another.
  int mantissa_bits;                        ///< Number of high bits to take from the mantissa in
                                            ///<  the table's corresponding floating point format
                                            ///<  when determining the table index
  ullint dp_detail_bitmask;                 ///< The logarithmic spline table divides a function
                                            ///<   argument, e.g. the squared distance between two
                                            ///<   particles, into an index (the high bits) and
                                            ///<   remaining detail (the low bits).  The nuances of
                                            ///<   the low bits must be translated by the spline's
                                            ///<   basis functions into an accurate estimate of the
                                            ///<   original function's true value.  This mask will
                                            ///<   isolate those bits from a double-precision real
                                            ///<   number.
  uint sp_detail_bitmask;                   ///< Single-precision variant of dp_detail_bitmask
  int segments_per_stride;                  ///< The number of spline segments sharing the same
                                            ///<   support width in one group for a single power of
                                            ///<   two within the overall range of the table.  This
                                            ///<   is 2^mantissa_bits.
  double ewald_coefficient;                 ///< The Ewald coefficient, when emulating PME direct
                                            ///<   space functions
  double coulomb_constant;                  ///< The coulomb constant to use in formulating any of
                                            ///<   the electrostatic PME interaction potentials
  double maximum_range;                     ///< The maximum range of the spline (this will affect
                                            ///<   the extent of the table)
  double minimum_absolute_range;            ///< The minimum absolute range for which splines will
                                            ///<   be computed.  This cannot be less than a range
                                            ///<   at which increments of the function argument
                                            ///<   computed in double precision would fall below
                                            ///<   machine precision, but choosing a larger value,
                                            ///<   i.e. 1.0e-5, at which the value of the tabulated
                                            ///<   function might be so high as to be irrelevant,
                                            ///<   is still feasible.  The value of the function at
                                            ///<   these minimal ranges will be assumed to remain
                                            ///<   constant beneath the minimum_absolute_range.
  double minimum_significant_range;         ///< The minimum significant range to consider when
                                            ///<   analyzing overall error in the spline.  The
                                            ///<   smallest ranges of the table become extremely
                                            ///<   dense and are likely to have the largest errors,
                                            ///<   but the chances of those indices being accessed
                                            ///<   are miniscule.  This is also the lower limit, in
                                            ///<   single-precision spline tables, at which ULP
                                            ///<   (units of least place) optimization commences.
  float indexing_offset;                    ///< Offset to be applied to whatever number is used to
                                            ///<   choose the index of the table.  This can clip
                                            ///<   off a vast number of the spline coefficient sets
                                            ///<   when the function is well behaved at short range
                                            ///<   and numbers near zero could be valid for
                                            ///<   indexing the table (which could cause cache
                                            ///<   thrashing).
  int ulp_optimization_depth;               ///< The depth of optimization to attempt in refining
                                            ///<   the units of least place.  Each floating-point
                                            ///<   number in the spline will be altered by up to
                                            ///<   2^ulp_optimization_depth increments of the value
                                            ///<   of its smallest bit.
  PrecisionModel precision;                 ///< Precision model for which the table is intended to
                                            ///<   operate (this is the "native" type, the
                                            ///<   corresponding type of the four-tuples from which
                                            ///<   the table is constructed)
  Hybrid<T4> table;                         ///< Array of spline coefficients for the logarithmic
                                            ///<   table.  There are 2^mantissa_bits splines for
                                            ///<   each range [ 2^N, 2^N+1 ) where N covers the
                                            ///<   range of the native type's exponent, less any
                                            ///<   upper extremes that the developer deems
                                            ///<   unimportant.
  double mean_overall_error;                ///< The average overall error found throughout the
                                            ///<   significant range of the spline table.
  double stdev_overall_error;               ///< Standard deviation of error found throughtout the
                                            ///<   spline table, over the range of interest
  double max_overall_error;                 ///< Maximum overall error identified throughout the
                                            ///<   spline, on the range of interest
  std::vector<double> mean_segment_error;   ///< The mean error estimated for individual spline
                                            ///<   segments (indices of this array track those of
                                            ///<   the table member variable)
  std::vector<double> stdev_segment_error;  ///< The standard deviations of error estimated for
                                            ///<   individual spline segments
  std::vector<double> max_segment_error;    ///< Maximum errors identified in each spline segment

  /// \brief Set the precision model.  Check the underlying data type in the process.
  void setPrecisionModel();
  
  /// \brief Set the target form of the function and check for the existence of the proper inputs
  ///
  /// \param target_form_in  The form of the function that the spline table is to emulate
  /// \param custom_form_in  Pre-computed tables with { argument, function value } pairs that can
  ///                        be queried to obtain the splines
  void setTargetForm(LogSplineForm target_form_in,
                     const std::vector<std::vector<double2>> &custom_form_in = {});

  /// \brief Set the Ewald coefficient after making validity checks.
  ///
  /// \param ewald_coefficient_in  The Ewald coefficient to set, equal to 1 / (2 sigma), where
  ///                              sigma is the Gaussian spread of charge or other density
  /// \param policy                Action to take in response to bad input
  void setEwaldCoefficient(double ewald_coefficient_in,
                           ExceptionResponse policy = ExceptionResponse::DIE);

  /// \brief Set the number of bits of the mantissa to use in indexing to the table.
  ///
  /// \param mantissa_bits_in  The number of bits in the mantissa to use in constructing the
  ///                          table index, beneath the exponent bits
  void setMantissaBits(int mantissa_bits_in,
                       ExceptionResponse policy = ExceptionResponse::WARN);

  /// \brief Set the minimum absolute range for which spline coefficients will be explicitly
  ///        computed.  The function imposes a lower bound based on the smallest increments that
  ///        an IEEE 64-bit floating point number can represent.  Beneath this range, spline
  ///        coefficients will copy the coefficients for the minimum interval which was explicitly
  ///        determined.
  ///
  /// \param min_range_in  The minimum range to apply
  void setMinimumAbsoluteRange(double min_range_in);

  /// \brief Set the number of optimization depth in terms of the units of least place (ULPs) for
  ///        each coefficient.
  ///
  /// \param ulp_optimization_depth_in  The portions of each cofficient's ULP to apply in trying to
  ///                                   optimize the spline results
  void setULPDepth(int ulp_optimization_depth_in);
  
  /// \brief Check the indexing offset to ensure its validity (that is is non-negative and a
  ///        power of two).  Return the input as the result if valid.
  ///
  /// \param indexing_offset_in  The requested indexing offset
  float checkIndexingOffset(float indexing_offset_in);
  
  /// \brief Set the size of the table, allocate, and modify the maximum range as appropriate.
  ///
  /// \brief max_range_in    The specified maximum range of the table
  /// \param custom_form_in  Pre-computed tables with { argument, function value } pairs that can
  ///                        be queried to obtain the splines
  /// \param policy          Action to take in response to bad input
  void allocate(double max_range_in, const std::vector<std::vector<double2>> &custom_form_in = {},
                const ExceptionResponse policy = ExceptionResponse::WARN);

  /// \brief Evaluate the spline table based on a constant rolled in.
  void evaluatePrescribedTable();

  /// \brief Compute the spline table due to a custom, pre-tabulated function.  This is a means for
  ///        reducing the size of a potential table, if it is well behaved at long range and very
  ///        steep at short range.
  void evaluateCustomTable();

  /// \brief Compute the function at a series of points over a segment.
  ///
  /// \param bwrk   Values of the function evaluated at even intervals, modified and returned
  /// \param bref   Reference values, a copy of bwrk.  If set to nullptr these values will not
  ///               be recorded.
  /// \param r2_i   The leftmost point of the segment, whether defined in terms of the squared
  ///               argument of the original function or the actual argument of the original
  ///               function
  /// \param width  The width of the segment in question
  void evaluateTargetInSegment(std::vector<double> *bwrk, std::vector<double> *bref, double r2_i,
                               double width);

  /// \brief Estimate the function at a series of points over a segment.
  ///
  /// \param estm     Values of the function evaluated at even intervals, modified and returned
  /// \param spl_tmp  Current guess for the spline parameter set
  /// \param r2_i     The leftmost point of the segment, whether defined in terms of the squared
  ///                 argument of the original function or the actual argument of the original
  ///                 function.  This is used when evaluating splines based on MIXED_FRACTION
  ///                 basis functions but not POLYNOMIAL basis functions.
  /// \param width    The width of the segment in question, used with evaluating splines based on
  ///                 MIXED_FRACTION basis functions
  void evaluateSplineInSegment(std::vector<double> *estm, const T4 spl_tmp, double r2_i,
                               double width);

  /// \brief Evaluate the errors over the range of splines deemed most relevant.
  void evaluateOverallError();
};

/// \brief Create the detail mask for a double-precision logarithmic spline table.
///
/// \param mantissa_bits  The number of bits, in addition to the exponent and sign bits, used to
///                       encode the table index
ullint doublePrecisionSplineDetailMask(int mantissa_bits);

/// \brief Create the detail mask for a single-precision logarithmic spline table.
///
/// \param mantissa_bits  The number of bits, in addition to the exponent and sign bits, used to
///                       encode the table index
uint singlePrecisionSplineDetailMask(int mantissa_bits);

/// \brief Compute the modified argument based on the "detail bits" for use in calculations with
///        POLYNOMIAL spline basis functions with DOUBLE-precision spline tables.
///
/// \param arg            The basic argument used to generate a spline table lookup index, e.g.
///                       the squared distance between two particles
/// \param mantissa_bits  The number of bits, in addition to the exponent and sign bits, used to
///                       encode the table index
/// \param detail_mask    The mask used to transform arg into the real-valued spline argument
double doublePrecisionSplineArgument(double arg, int mantissa_bits, ullint detail_mask);

/// \brief Compute the modified argument based on the "detail bits" for use in calculations with
///        POLYNOMIAL spline basis functions with SINGLE-precision spline tables.  Descriptions of
///        input arguments follow from doublePrecisionSplineArgument(), above.
float singlePrecisionSplineArgument(float arg, int mantissa_bits, uint detail_mask);

/// \brief Restore the templated type of a LogScaleSpline object's abstract.  Overloads of this
///        function handle type restoration in their respective objects.
///
/// \param rasa  The "blank slate" LogScaleSpline abstract, with no templated type characteristics
/// \{
template <typename T4> LogSplineTable<T4> restoreType(const LogSplineTable<void> *rasa);

template <typename T4> LogSplineTable<T4> restoreType(const LogSplineTable<void> &rasa);
/// \}
  
} // namespace stmath
} // namespace stormm

#include "log_scale_spline.tpp"

#endif
