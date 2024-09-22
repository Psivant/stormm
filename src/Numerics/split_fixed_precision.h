// -*-c++-*-
#ifndef STORMM_SPLIT_FIXED_PRECISION_H
#define STORMM_SPLIT_FIXED_PRECISION_H

#include <string>
#include <climits>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "Constants/scaling.h"
#include "DataTypes/stormm_vector_types.h"
#include "numeric_enumerators.h"

namespace stormm {
namespace numerics {

using card::Hybrid;
using constants::ExceptionResponse;
using constants::PrecisionModel;
  
/// \brief The maximum contributions for signed integer accumulation.  There is no long long
///        integer form of the maximum long long integer accumulation as the number would wrap
///        the format to become the negative of its intended value.  The long long integer
///        representation of the maximum int accumulation is present to facilitate conversion to
///        a unified 64-bit integer, but even that is only reliable on the CPU (the NVIDIA CUDA
///        compiler seems to take type specifications of long long int as a suggestion, not a
///        command, and this may relate to thread counts and register pressure as the compiler
///        tries to optimize a kernel for given launch bounds).
/// \{
constexpr llint max_int_accumulation_ll = (1LL << (int_bit_count_int - 1));
constexpr double max_int_accumulation   = max_int_accumulation_ll;
constexpr float max_int_accumulation_f  = max_int_accumulation;
constexpr double max_llint_accumulation = max_int_accumulation * max_int_accumulation * 2.0;
constexpr float max_llint_accumulation_f  = max_llint_accumulation;
constexpr int max_short_accumulation = (1 << 15);
/// \}

/// \brief Translate a string specifying a force accumulation method into the numerical code.
///
/// \param method  The string to translate
/// \param policy  Action to take in the event of an error
AccumulationMethod translateAccumulationMethod(const std::string &choice,
                                               ExceptionResponse policy);

/// \brief Get a string for the name of a force accumulation method.
///
/// \param method  The method in question
std::string getAccumulationMethodName(AccumulationMethod method);

/// \brief Determine the best accumulation method based on the precision level of the forces.
///
/// \param frc_bits  Number of bits stored after the decimal in fixed-precision force
///                  representations
AccumulationMethod chooseAccumulationMethod(int frc_bits);

/// \brief Produce an error message describing range violations in user choices for various
///        fixed-precision methods.
///
/// \param choice   The selected number of bits in the precision model
/// \param min_val  The minimum allowed number of bits
/// \param max_val  The maximum allowed number of bits
std::string fixedPrecisionRangeErrorMessage(int choice, int min_val, int max_val);

/// \brief Check user input regarding the global position scaling.  Positions are represented in
///        internal units of Angstroms.
///
/// \param choice  The fixed-precision bits for representing global positions of particles
void checkGlobalPositionBits(int choice);

/// \brief Check user input regarding the local position scaling.  Positions are represented in
///        internal units of Angstroms.
///
/// \param choice  The fixed-precision bits for representing local positions of particles
void checkLocalPositionBits(int choice);

/// \brief Check user input regarding the fixed-precision velocity scaling.  Velocities are
///        represented in units of Angstroms per sqrt(418.4) * femtoseconds.
///
/// \param choice  The fixed-precision bits for representing particle velocities
void checkVelocityBits(int choice);

/// \brief Check user input regarding the fixed-precision force accumulation.  Forces are
///        represented in units of kcal/mol - Angstroms.
///
/// \param choice  The fixed-precision bits for representing forces acting on particles
void checkForceBits(int choice);

/// \brief Check user input regarding the fixed-precision energy accumulation.  Energies are
///        represented in units of kcal/mol.
///
/// \param choice  The fixed-precision bits for representing energy contributions
void checkEnergyBits(int choice);

/// \brief Check user input regarding the fixed-precision charge denisty accumulation on the mesh.
///        Charge density per grid point is represented in atomic units.
///
/// \param choice  The fixed-precision bits for performing charge density accumulation
/// \param pmdoel  The fixed-precision model, which implies the accumulation range
void checkChargeMeshBits(int choice, PrecisionModel pmodel);

/// \brief Consolidate repetitive vector size comparisons and errors messages in the following
///        conversion functions.
///
/// \param primary_length   Length of the primary vector
/// \param overflow_length  Length of the overflow vector
/// \param caller           Name of the calling function
/// \param message          The error message to print
void checkFPVectorLength(size_t primary_length, size_t overflow_length, const char* caller,
                         const std::string &message = std::string(""));

/// \brief Convert floating point numbers into fixed-precision representations with two integers.
///        This is similar to splitRealAccumulation below, but will set the values rather than
///        add new contributions.
///
/// Overloaded:
///   - Convert a single-precision floating point number into two 32-bit signed integers.
///   - Convert a double-precision floating point number into a 64-bit primary integer and a
///     32-bit secondary / overflow integer.
///   - Convert N values in a single array to one pair of integer arrays
///   - Convert N values in three arays to three pairs of integer arrays
///   - Return the fixed-precision representation or assign it directly to the corresponding
///     locations in two appropriate arrays.
///   - When working with arrays, accept C-style arrays, Standard Template Library Vectors, or
///     Hybrid objects (the conversion must be done on the HOST side)
///
/// \param fval      Single-precision value to convert to fixed-precision
/// \param dval      Double-precision value to convert to fixed-precision
/// \param primary   The primary accumulator (the low 32 bits)
/// \param overflow  The secondary accumulator (the high 31 bits)
/// \param n_values  The number of values in any of the original real-number arrays
/// \param scale     When submitting arrays for conversion, this
/// \{
int2 hostFloatToInt63(const float fval);

void hostFloatToInt63(const float fval, int *primary, int *overflow);

void hostFloatToInt63(const float* fval, int* primary, int* overflow, size_t n_values,
                      float scale = 1.0f);

void hostFloatToInt63(const float* fval_x, const float* fval_y, const float* fval_z,
                      int* primary_x, int* overflow_x, int* primary_y, int* overflow_y,
                      int* primary_z, int* overflow_z, size_t n_values, float scale = 1.0f);

void hostFloatToInt63(const std::vector<float> &fval, std::vector<int> *primary,
                      std::vector<int> *overflow, float scale = 1.0f);

void hostFloatToInt63(const Hybrid<float> &fval, Hybrid<int> *primary, Hybrid<int> *overflow,
                      float scale = 1.0f);

void hostFloatToInt63(const std::vector<float> &fval_x, const std::vector<float> &fval_y,
                      const std::vector<float> &fval_z, std::vector<int> *primary_x,
                      std::vector<int> *overflow_x, std::vector<int> *primary_y,
                      std::vector<int> *overflow_y, std::vector<int> *primary_z,
                      std::vector<int> *overflow_z, float scale = 1.0f);

void hostFloatToInt63(const Hybrid<float> &fval_x, const Hybrid<float> &fval_y,
                      const Hybrid<float> &fval_z, Hybrid<int> *primary_x, Hybrid<int> *overflow_x,
                      Hybrid<int> *primary_y, Hybrid<int> *overflow_y, Hybrid<int> *primary_z,
                      Hybrid<int> *overflow_z, float scale = 1.0f);

int2 hostLongLongToInt63(const llint val);

void hostLongLongToInt63(const llint val, int *primary, int *overflow);

void hostLongLongToInt63(const llint* val, int* primary, int* overflow, size_t n_values);

void hostLongLongToInt63(const llint* val_x, const llint* val_y, const llint* val_z,
                         int* primary_x, int* overflow_x, int* primary_y, int* overflow_y,
                         int* primary_z, int* overflow_z, size_t n_values);

void hostLongLongToInt63(const std::vector<llint> &val, std::vector<int> *primary,
                         std::vector<int> *overflow);

void hostLongLongToInt63(const std::vector<llint> &val_x, const std::vector<llint> &val_y,
                         const std::vector<llint> &val_z, std::vector<int> *primary_x,
                         std::vector<int> *overflow_x, std::vector<int> *primary_y,
                         std::vector<int> *overflow_y, std::vector<int> *primary_z,
                         std::vector<int> *overflow_z);

void hostLongLongToInt63(const Hybrid<llint> &val, Hybrid<int> *primary, Hybrid<int> *overflow);

void hostLongLongToInt63(const Hybrid<llint> &val_x, const Hybrid<llint> &val_y,
                         const Hybrid<llint> &val_z, Hybrid<int> *primary_x,
                         Hybrid<int> *overflow_x, Hybrid<int> *primary_y, Hybrid<int> *overflow_y,
                         Hybrid<int> *primary_z, Hybrid<int> *overflow_z);

int2 hostDoubleToInt63(const double fval);

void hostDoubleToInt63(const double fval, int *primary, int *overflow);

int95_t hostDoubleToInt95(const double fval);

void hostDoubleToInt95(const double fval, llint *primary, int *overflow);

void hostDoubleToInt95(const double* dval, llint* primary, int* overflow, size_t n_values,
                       double scale = 1.0);

void hostDoubleToInt95(const std::vector<double> &dval, std::vector<llint> *primary,
                       std::vector<int> *overflow, double scale = 1.0);

void hostDoubleToInt95(const Hybrid<double> &dval, Hybrid<llint> *primary, Hybrid<int> *overflow,
                       double scale = 1.0);
  
void hostDoubleToInt95(const double* dval_x, const double* dval_y, const double* dval_z,
                       llint* primary_x, int* overflow_x, llint* primary_y, int* overflow_y,
                       llint* primary_z, int* overflow_z, size_t n_values, double scale = 1.0);

void hostDoubleToInt95(const std::vector<double> &dval_x, const std::vector<double> &dval_y,
                       const std::vector<double> &dval_z, std::vector<llint> *primary_x,
                       std::vector<int> *overflow_x, std::vector<llint> *primary_y,
                       std::vector<int> *overflow_y, std::vector<llint> *primary_z,
                       std::vector<int> *overflow_z, double scale = 1.0);

void hostDoubleToInt95(const Hybrid<double> &dval_x, const Hybrid<double> &dval_y,
                       const Hybrid<double> &dval_z, Hybrid<llint> *primary_x,
                       Hybrid<int> *overflow_x, Hybrid<llint> *primary_y, Hybrid<int> *overflow_y,
                       Hybrid<llint> *primary_z, Hybrid<int> *overflow_z, double scale = 1.0);
/// \}

/// \brief Convert numbers in split fixed precision to floating point reals.  Downscaling to the
///        proper units is the responsibility of the developer.
///
/// Overloaded:
///   - Convert two 32-bit integer values or a 64-bit / 32-bit combination into a float or double
///     (64-bit primary accumulators will always convert to double, but the output can be recast
///     as float)
///   - Convert the corresponding fused tuples into either type (int95_t will be converted to
///     double only)
///   - Convert one array of numbers, or even three arrays, of a trusted length
///   - Convert split fixed-precision integers in C-style arrays, Standard Template Library
///     vectors, or Hybrid objects
///
/// \param primary   Primary accumulator (the maximum accumulation increment is inferred from the
///                  size of the data type)
/// \param overflow  Overflow accumulator
/// \param result    Array to collect the results of the conversion
/// \param n_values  Trusted length of result, as well as any primary and overflow arrays, when
///                  working with C-style arrays in the conversion
/// \{
llint hostInt63ToLongLong(int primary, int overflow);

void hostInt63ToLongLong(llint* result, const int* primary, const int* overflow, size_t n_values);

void hostInt63ToLongLong(std::vector<llint> *result, const std::vector<int> &primary,
                         const std::vector<int> &overflow, size_t n_values);

void hostInt63ToLongLong(Hybrid<llint> *result, const Hybrid<int> &primary,
                         const Hybrid<int> &overflow, size_t n_values);

void hostInt63ToLongLong(llint* result_x, llint* result_y, llint* result_z, const int* primary_x,
                         const int* overflow_x, const int* primary_y, const int* overflow_y,
                         const int* primary_z, const int* overflow_z, size_t n_values);

void hostInt63ToLongLong(std::vector<llint> *result_x, std::vector<llint> *result_y,
                         std::vector<llint> *result_z, const std::vector<int> &primary_x,
                         const std::vector<int> &overflow_x, const std::vector<int> &primary_y,
                         const std::vector<int> &overflow_y, const std::vector<int> &primary_z,
                         const std::vector<int> &overflow_z, size_t n_values);

void hostInt63ToLongLong(Hybrid<llint> *result_x, Hybrid<llint> *result_y, Hybrid<llint> *result_z,
                         const Hybrid<int> &primary_x, const Hybrid<int> &overflow_x,
                         const Hybrid<int> &primary_y, const Hybrid<int> &overflow_y,
                         const Hybrid<int> &primary_z, const Hybrid<int> &overflow_z,
                         size_t n_values);

double hostInt63ToDouble(int primary, int overflow);

void hostInt63ToDouble(double* result, const int* primary, const int* overflow, size_t n_values,
                       double descale = 1.0);

void hostInt63ToDouble(std::vector<double> *result, const std::vector<int> &primary,
                       const std::vector<int> &overflow, double descale = 1.0);

void hostInt63ToDouble(Hybrid<double> *result, const Hybrid<int> &primary,
                       const Hybrid<int> &overflow, double descale = 1.0);

float hostInt63ToFloat(int primary, int overflow);

void hostInt63ToFloat(float* result, const int* primary, const int* overflow, size_t n_values,
                      float descale = 1.0f);

void hostInt63ToFloat(std::vector<float> *result, const std::vector<int> &primary,
                      const std::vector<int> &overflow, float descale = 1.0f);

void hostInt63ToFloat(Hybrid<float> *result, const Hybrid<int> &primary,
                      const Hybrid<int> &overflow, float descale = 1.0f);

void hostInt63ToDouble(double* result_x, double* result_y, double* result_z, const int* primary_x,
                       const int* overflow_x, const int* primary_y, const int* overflow_y,
                       const int* primary_z, const int* overflow_z, size_t n_values,
                       double descale = 1.0);

void hostInt63ToDouble(std::vector<double> *result_x, std::vector<double> *result_y,
                       std::vector<double> *result_z, const std::vector<int> &primary_x,
                       const std::vector<int> &overflow_x, const std::vector<int> &primary_y,
                       const std::vector<int> &overflow_y, const std::vector<int> &primary_z,
                       const std::vector<int> &overflow_z, size_t n_values, double descale = 1.0);

void hostInt63ToDouble(Hybrid<double> *result_x, Hybrid<double> *result_y,
                       Hybrid<double> *result_z, const Hybrid<int> &primary_x,
                       const Hybrid<int> &overflow_x, const Hybrid<int> &primary_y,
                       const Hybrid<int> &overflow_y, const Hybrid<int> &primary_z,
                       const Hybrid<int> &overflow_z, size_t n_values, double descale = 1.0);

void hostInt63ToFloat(float* result_x, float* result_y, float* result_z, const int* primary_x,
                      const int* overflow_x, const int* primary_y, const int* overflow_y,
                      const int* primary_z, const int* overflow_z, size_t n_values,
                      float descale = 1.0f);

void hostInt63ToFloat(std::vector<float> *result_x, std::vector<float> *result_y,
                      std::vector<float> *result_z, const std::vector<int> &primary_x,
                      const std::vector<int> &overflow_x, const std::vector<int> &primary_y,
                      const std::vector<int> &overflow_y, const std::vector<int> &primary_z,
                      const std::vector<int> &overflow_z, float descale = 1.0f);

void hostInt63ToFloat(Hybrid<float> *result_x, Hybrid<float> *result_y, Hybrid<float> *result_z,
                      const Hybrid<int> &primary_x, const Hybrid<int> &overflow_x,
                      const Hybrid<int> &primary_y, const Hybrid<int> &overflow_y,
                      const Hybrid<int> &primary_z, const Hybrid<int> &overflow_z,
                      float descale = 1.0f);

double hostInt63ToDouble(int2 ival);

float hostInt63ToFloat(int2 ival);

double hostInt95ToDouble(const int95_t ival);
  
double hostInt95ToDouble(llint primary, int overflow);

void hostInt95ToDouble(double* result, const llint* primary, const int* overflow, size_t n_values,
                       double descale = 1.0);

void hostInt95ToDouble(std::vector<double> *result, const std::vector<llint> &primary,
                       const std::vector<int> &overflow, double descale = 1.0);

void hostInt95ToDouble(Hybrid<double> *result, const Hybrid<llint> &primary,
                       const Hybrid<int> &overflow, double descale = 1.0);

void hostInt95ToDouble(double* result_x, double* result_y, double* result_z,
                       const llint* primary_x, const int* overflow_x, const llint* primary_y,
                       const int* overflow_y, const llint* primary_z, const int* overflow_z,
                       size_t n_values, double descale = 1.0);

void hostInt95ToDouble(std::vector<double> *result_x, std::vector<double> *result_y,
                       std::vector<double> *result_z, const std::vector<llint> &primary_x,
                       const std::vector<int> &overflow_x, const std::vector<llint> &primary_y,
                       const std::vector<int> &overflow_y, const std::vector<llint> &primary_z,
                       const std::vector<int> &overflow_z, double descale = 1.0);

void hostInt95ToDouble(Hybrid<double> *result_x, Hybrid<double> *result_y,
                       Hybrid<double> *result_z, const Hybrid<llint> &primary_x,
                       const Hybrid<int> &overflow_x, const Hybrid<llint> &primary_y,
                       const Hybrid<int> &overflow_y, const Hybrid<llint> &primary_z,
                       const Hybrid<int> &overflow_z, double descale = 1.0);
/// \}

/// \brief Accumulate floating point numbers into fixed-precision representations with two
///        integers.
///
/// Overloaded:
///   - Convert a single-precision floating point number into two 32-bit signed integers.
///   - Convert a double-precision floating point number into a 64-bit primary integer and a
///     32-bit secondary / overflow integer.
///
/// \param fval      Single-precision value to convert to fixed-precision
/// \param dval      Double-precision value to convert to fixed-precision
/// \param primary   The primary accumulator (the low 32 bits)
/// \param overflow  The secondary accumulator (the high 31 bits)
/// \{
void hostSplitAccumulation(const float fval, int *primary, int *overflow);

void hostSplitAccumulation(const double fval, llint *primary, int *overflow);
/// \}

/// \brief Accumulate two split fixed-precision integers.  This routine is unsafe for summing one
///        split fixed-precision number with the "negative" of another, because { -a, -b } is not
///        -{ a, b } (where "a" is the primary accumulator and b the secondary accumulator) if "a"
///        is the minimum value of its respective integer, e.g. INT_MIN.  This is because -INT_MIN
///        would overflow, having a value one greater than the corresponding INT_MAX, equivalent
///        to stating -INT_MIN = INT_MIN.  To use this function (or equivalent HPC code in
///        accumulate.cui) to perform 'subtraction' will work in all other cases, but because the
///        effect of the overflow is so rare and can have effects ranging from subtle to
///        catastrophic such an approach must not be used.  The corresponding host(...)Subtract()
///        functions below handle the subtraction of one split fixed-precision number from another.
///
/// Overloaded:
///   - Accumulate int95_t with its various components
///   - Accumulate int63_t (int2) with its various components
///
/// \param a      The first of the two value pairs
/// \param b      The second of the two value pairs
/// \param a_x    The lower bits of the first value pair
/// \param a_y    The upper bits of the first value pair
/// \param b_x    The lower_bits of the second value pair
/// \param b_y    The upper_bits of the second value pair
/// \param breal  Real-valued form of the second number (this will be converted to the appropriate
///               split fixed-precision type before adding)
/// \{
int95_t hostSplitFPSum(const int95_t a, const int95_t b);

int2 hostSplitFPSum(const int2 a, const int2 b);

int95_t hostSplitFPSum(const int95_t a, double breal);

int2 hostSplitFPSum(const int2 a, float breal);

int95_t hostSplitFPSum(const int95_t a, llint b_x, int b_y);

int2 hostSplitFPSum(const int2 a, int b_x, int b_y);

int95_t hostInt95Sum(llint a_x, int a_y, llint b_x, int b_y);

int2 hostInt63Sum(int a_x, int a_y, int b_x, int b_y);

int95_t hostInt95Sum(llint a_x, int a_y, double breal);

int2 hostInt63Sum(int a_x, int a_y, float breal);
/// \}

/// \brief Subtract one split fixed-precision integer from another of the same type.  This includes
///        a guard against the case where the integer to be subtracted holds the minimum value,
///        e.g. INT_MIN, in its primary accumulator.  Overloading and descriptions of input
///        parameters follow from the equivalent host(...)Sum() functions above.
/// \{
int95_t hostSplitFPSubtract(const int95_t a, const int95_t b);

int2 hostSplitFPSubtract(const int2 a, const int2 b);

int95_t hostSplitFPSubtract(const int95_t a, double breal);

int2 hostSplitFPSubtract(const int2 a, float breal);

int95_t hostSplitFPSubtract(const int95_t a, llint b_x, int b_y);

int2 hostSplitFPSubtract(const int2 a, int b_x, int b_y);

int95_t hostInt95Subtract(llint a_x, int a_y, llint b_x, int b_y);

int2 hostInt63Subtract(int a_x, int a_y, int b_x, int b_y);

int95_t hostInt95Subtract(llint a_x, int a_y, double breal);

int2 hostInt63Subtract(int a_x, int a_y, float breal);
/// \}

/// \brief Multiply a split fixed-precision number by a 32-bit integer.
///
/// Overloaded:
///   - Multiplty a 63-bit split fixed-precision integer
///   - Multiplty a 95-bit split fixed-precision integer
///
/// \param a  The fused, split fixed-precision integer
/// \param b  The multiplying factor
/// \{
int95_t hostSplitFPMult(const int95_t a, int b);

int95_t hostInt95Mult(llint a_x, int a_y, int b);

int2 hostSplitFPMult(const int2 a, int b);

int2 hostInt63Mult(int a_x, int a_y, int b);
/// \}
  
/// \brief Convert a split fixed-precision number with one bit scaling into an equivalent data type
///        with a different bit scaling.
///
/// Overloaded:
///   - Convert int63_t (int2) to another int63_t
///   - Convert int95_t to another int95_t
///   - Convert a Standard Template Library vector of int95_t
///   - Convert two C-style arrays with a trusted length, two Standard Template Library vectors, or
///     two Hybrid objects of the appropriate data types
///
/// \param fp           The fixed-precision number of interest, or the primary array of values
/// \param fp_ovrf      Overflow bits to complement fp
/// \param length       Trusted length of C-style array inputs
/// \param native_bits  Original bit scaling of the fixed-precision format
/// \param output_bits  Bit scaling of the output format
/// \{
int2 hostChangeFPBits(const int2 fp, int native_bits, int output_bits);

int95_t hostChangeFPBits(const int95_t fp, int native_bits, int output_bits);

void hostChangeFPBits(std::vector<int95_t> *fp, int native_bits, int output_bits);

void hostChangeFPBits(int* fp, int* fp_ovrf, size_t length, int native_bits, int output_bits);

void hostChangeFPBits(llint* fp, int* fp_ovrf, size_t length, int native_bits, int output_bits);

void hostChangeFPBits(std::vector<int> *fp, std::vector<int> *fp_ovrf, int native_bits,
                      int output_bits);

void hostChangeFPBits(std::vector<llint> *fp, std::vector<int> *fp_ovrf, int native_bits,
                      int output_bits);

void hostChangeFPBits(Hybrid<int> *fp, Hybrid<int> *fp_ovrf, int native_bits, int output_bits);

void hostChangeFPBits(Hybrid<llint> *fp, Hybrid<int> *fp_ovrf, int native_bits, int output_bits);
/// \}

/// \brief Compute a one-dimensional grid of points, beginning at some arbitrary origin and
///        marking at uniform increments as the grid line evolves in the Cartesian X, Y, and Z
///        dimensions.
///
/// Overloaded:
///   - Populate an Standard Template Library vectorof int95_t data
///   - Populate an llint (primary) vector and an int (overflow) vector
///   - Populate the HOST data of llint (primary) and int (overflow) Hybrid objects
///   - Populate llint (primary) and int (overflow) C-style arrays
///
/// \param coordinates  The array on which to accumulate the grid.  Filled and returned.
/// \param primary      Primary bits for the fixed-precision accumulators.  Filled and returned.
/// \param overflow     Overflow bits for the fixed-precision accumulators.  Filled and returned.
/// \param origin       The origin of the mesh and location of is first point
/// \param increment    The spacing between grid points
/// \param grid_size    Provided for C-style arrays, defining their trusted length
/// \{
void fixedPrecisionGrid(std::vector<int95_t> *coordinates, const int95_t origin,
                        const int95_t increment);

void fixedPrecisionGrid(std::vector<llint> *primary, std::vector<int> *overflow,
                        const int95_t origin, const int95_t increment);

void fixedPrecisionGrid(Hybrid<llint> *primary, Hybrid<int> *overflow, const int95_t origin,
                        const int95_t increment);

void fixedPrecisionGrid(llint *primary, int *overflow, const int95_t origin,
                        const int95_t increment, const size_t grid_size);
/// \}
  
} // namespace numerics
} // namespace stormm

// Make basic arithmetic functions native to the stormm namespace
namespace stormm {
  using numerics::hostFloatToInt63;
  using numerics::hostDoubleToInt95;
  using numerics::hostInt63ToDouble;
  using numerics::hostInt63ToFloat;
  using numerics::hostInt95ToDouble;
  using numerics::hostSplitAccumulation;
  using numerics::hostSplitFPSubtract;
  using numerics::hostSplitFPSum;
  using numerics::hostInt95Subtract;
  using numerics::hostInt95Sum;
  using numerics::hostInt63Subtract;
  using numerics::hostInt63Sum;
  using numerics::max_int_accumulation;
  using numerics::max_int_accumulation_f;
  using numerics::max_int_accumulation_ll;
  using numerics::max_llint_accumulation;
  using numerics::max_llint_accumulation_f;
  using numerics::max_short_accumulation;
} // namespace stormm

#endif

