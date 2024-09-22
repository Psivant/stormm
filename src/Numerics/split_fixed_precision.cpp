#include "copyright.h"
#include "Constants/fixed_precision.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "Reporting/error_format.h"
#include "split_fixed_precision.h"

namespace stormm {
namespace numerics {

using constants::getEnumerationName;
using parse::strcmpCased;
using parse::PolyNumeric;

//-------------------------------------------------------------------------------------------------
AccumulationMethod translateAccumulationMethod(const std::string &choice,
                                                         const ExceptionResponse policy) {
  if (strcmpCased(choice, std::string("split"))) {
    return AccumulationMethod::SPLIT;
  }
  else if (strcmpCased(choice, std::string("whole"))) {    
    return AccumulationMethod::WHOLE;
  }
  else if (strcmpCased(choice, std::string("automatic"))) {    
    return AccumulationMethod::AUTOMATIC;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getAccumulationMethodName(const AccumulationMethod method) {
  switch (method) {
  case AccumulationMethod::SPLIT:
    return std::string("SPLIT");
  case AccumulationMethod::WHOLE:
    return std::string("WHOLE");
  case AccumulationMethod::AUTOMATIC:
    return std::string("AUTOMATIC");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
AccumulationMethod chooseAccumulationMethod(const int frc_bits) {
  if (frc_bits <= 24) {
    return AccumulationMethod::SPLIT;
  }
  else {
    return AccumulationMethod::WHOLE;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string fixedPrecisionRangeErrorMessage(const int choice, const int min_val,
                                            const int max_val) {
  if (choice < min_val) {
    return std::string(" is too small, and would result in an inaccurate simulation.");
  }
  if (choice > max_val) {
    return std::string(" is too large, and might overflow the 64-bit integer format.");
  }
  return std::string("");
}

//-------------------------------------------------------------------------------------------------
void checkGlobalPositionBits(const int choice) {
  if (choice < min_globalpos_scale_bits || choice > max_globalpos_scale_bits) {
    rtErr("Global position scaling is allowed in a fixed precision range of " +
          std::to_string(min_globalpos_scale_bits) + " to " +
          std::to_string(max_globalpos_scale_bits) + ".  A value of " +
          std::to_string(choice) +
          fixedPrecisionRangeErrorMessage(choice, min_globalpos_scale_bits,
                                          max_globalpos_scale_bits), "checkGlobalPositionBits");
  }
}

//-------------------------------------------------------------------------------------------------
void checkLocalPositionBits(const int choice) {
  if (choice < min_localpos_scale_bits || choice > max_localpos_scale_bits) {
    rtErr("Local position scaling is allowed in a fixed precision range of " +
          std::to_string(min_localpos_scale_bits) + " to " +
          std::to_string(max_localpos_scale_bits) + " bits.  A value of " +
          std::to_string(choice) + fixedPrecisionRangeErrorMessage(choice, min_localpos_scale_bits,
                                                                   max_localpos_scale_bits),
          "checkLocalPositionBits");
  }
}

//-------------------------------------------------------------------------------------------------
void checkVelocityBits(const int choice) {
  if (choice < min_velocity_scale_bits || choice > max_velocity_scale_bits) {
    rtErr("Velocity scaling is allowed in a fixed precision range of " +
          std::to_string(min_velocity_scale_bits) + " to " +
          std::to_string(max_velocity_scale_bits) + " bits.  A value of " +
          std::to_string(choice) + fixedPrecisionRangeErrorMessage(choice, min_velocity_scale_bits,
                                                                   max_velocity_scale_bits),
          "checkVelocityBits");
  }
}

//-------------------------------------------------------------------------------------------------
void checkForceBits(const int choice) {
  if (choice < min_force_scale_bits || choice > max_force_scale_bits) {
    rtErr("Force accumulation is allowed in a fixed precision range of " +
          std::to_string(min_force_scale_bits) + " to " +
          std::to_string(max_force_scale_bits) + " bits.  A value of " + std::to_string(choice) +
          fixedPrecisionRangeErrorMessage(choice, min_force_scale_bits, max_force_scale_bits),
          "checkForceBits");
  }
}

//-------------------------------------------------------------------------------------------------
void checkEnergyBits(const int choice) {
  if (choice < min_energy_scale_bits || choice > max_energy_scale_bits) {
    rtErr("Energy accumulation is allowed in a fixed precision range of " +
          std::to_string(min_energy_scale_bits) + " to " +
          std::to_string(max_energy_scale_bits) + " bits.  A value of " + std::to_string(choice) +
          fixedPrecisionRangeErrorMessage(choice, min_energy_scale_bits, max_energy_scale_bits),
          "checkEnergyBits");
  }
}

//-------------------------------------------------------------------------------------------------
void checkChargeMeshBits(const int choice, const PrecisionModel pmodel) {
  if (choice < min_charge_mesh_scale_bits) {
    rtErr("Charge mesh accumulation must take place with at least " +
          std::to_string(min_charge_mesh_scale_bits) + " bits.  A values of " +
          std::to_string(choice) +
          fixedPrecisionRangeErrorMessage(choice, min_charge_mesh_scale_bits, 64),
          "checkChargeMeshBits");
  }
  switch (pmodel) {
  case PrecisionModel::SINGLE:
    if (choice > 31) {
      rtErr("Charge mesh accumulation in a " + getEnumerationName(pmodel) + " precision model "
            "will take place with a 32-bit signed integer accumulation grid.  A value of " +
            std::to_string(choice) + fixedPrecisionRangeErrorMessage(choice, 8, 31),
            "checkChargeMeshBits");
    }
    break;
  case PrecisionModel::DOUBLE:
    if (choice > max_charge_mesh_scale_bits) {
      rtErr("Charge mesh accumulation in a " + getEnumerationName(pmodel) + " precision model "
            "will take place with a 64-bit signed integer accumulation grid.  A value of " +
            std::to_string(choice) +
            fixedPrecisionRangeErrorMessage(choice, 8, max_charge_mesh_scale_bits),
            "checkChargeMeshBits");
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void checkFPVectorLength(const size_t primary_length, const size_t overflow_length,
                         const char* caller, const std::string &message) {
  if (primary_length != overflow_length) {
    std::string out_message;
    if (message.size() > 0) {
      out_message = message;
    }
    else {
      out_message = std::string("The primary and overflow vectors must be the same lengths (") +
                    std::to_string(primary_length) + " != " + std::to_string(overflow_length) +
                    ").";
    }
    rtErr(out_message, caller);
  }
}

//-------------------------------------------------------------------------------------------------
int2 hostFloatToInt63(const float fval) {
  int2 result;
  if (fabsf(fval) >= max_int_accumulation_f) {
    const int spillover = fval / max_int_accumulation_f;
    result.x = fval - (static_cast<float>(spillover) * max_int_accumulation_f);
    result.y = spillover;
  }
  else {
    result.x = fval;
    result.y = 0;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void hostFloatToInt63(const float fval, int *primary, int *overflow) {
  const int2 result = hostFloatToInt63(fval);
  *primary  = result.x;
  *overflow = result.y;
}

//-------------------------------------------------------------------------------------------------
void hostFloatToInt63(const float* fval, int* primary, int* overflow, const size_t n_values,
                      const float scale) {
  for (size_t i = 0LLU; i < n_values; i++) {
    hostFloatToInt63(fval[i] * scale, &primary[i], &overflow[i]);
  }
}

//-------------------------------------------------------------------------------------------------
void hostFloatToInt63(const std::vector<float> &fval, std::vector<int> *primary,
                      std::vector<int> *overflow, const float scale) {
  checkFPVectorLength(primary->size(), overflow->size(), "hostFloatToInt63");
  checkFPVectorLength(primary->size(), fval.size(), "hostFloatToInt63", "The original data must "
                      "have the same length as the integral result vectors.");
  hostFloatToInt63(fval.data(), primary->data(), overflow->data(), fval.size(), scale);
}

//-------------------------------------------------------------------------------------------------
void hostFloatToInt63(const Hybrid<float> &fval, Hybrid<int> *primary, Hybrid<int> *overflow,
                      const float scale) {
  checkFPVectorLength(primary->size(), overflow->size(), "hostFloatToInt63");
  checkFPVectorLength(primary->size(), fval.size(), "hostFloatToInt63", "The original data must "
                      "have the same length as the integral result vectors.");
  hostFloatToInt63(fval.data(), primary->data(), overflow->data(), fval.size(), scale);
}

//-------------------------------------------------------------------------------------------------
void hostFloatToInt63(const float* fval_x, const float* fval_y, const float* fval_z,
                      int* primary_x, int* overflow_x, int* primary_y, int* overflow_y,
                      int* primary_z, int* overflow_z, const size_t n_values, const float scale) {
  for (size_t i = 0LLU; i < n_values; i++) {
    hostFloatToInt63(fval_x[i] * scale, &primary_x[i], &overflow_x[i]);
    hostFloatToInt63(fval_y[i] * scale, &primary_y[i], &overflow_y[i]);
    hostFloatToInt63(fval_z[i] * scale, &primary_z[i], &overflow_z[i]);
  }
}

//-------------------------------------------------------------------------------------------------
void hostFloatToInt63(const std::vector<float> &fval_x, const std::vector<float> &fval_y,
                      const std::vector<float> &fval_z, std::vector<int> *primary_x,
                      std::vector<int> *overflow_x, std::vector<int> *primary_y,
                      std::vector<int> *overflow_y, std::vector<int> *primary_z,
                      std::vector<int> *overflow_z, const float scale) {
  hostFloatToInt63(fval_x, primary_x, overflow_x, scale);
  hostFloatToInt63(fval_y, primary_y, overflow_y, scale);
  hostFloatToInt63(fval_z, primary_z, overflow_z, scale);
}

//-------------------------------------------------------------------------------------------------
void hostFloatToInt63(const Hybrid<float> &fval_x, const Hybrid<float> &fval_y,
                      const Hybrid<float> &fval_z, Hybrid<int> *primary_x, Hybrid<int> *overflow_x,
                      Hybrid<int> *primary_y, Hybrid<int> *overflow_y, Hybrid<int> *primary_z,
                      Hybrid<int> *overflow_z, const float scale) {
  hostFloatToInt63(fval_x, primary_x, overflow_x, scale);
  hostFloatToInt63(fval_y, primary_y, overflow_y, scale);
  hostFloatToInt63(fval_z, primary_z, overflow_z, scale);
}

//-------------------------------------------------------------------------------------------------
int2 hostLongLongToInt63(const llint val) {
  int2 result;
  const int mult = (val < 0) ? -1 : 1;
  llint tval = val * mult;
  result.y = (tval >> 31);
  llint remainder = result.y;
  remainder <<= 31;
  remainder = tval - remainder;
  result.x = remainder * mult;
  result.y *= mult;
  return result;
}

//-------------------------------------------------------------------------------------------------
void hostLongLongToInt63(const llint val, int *primary, int *overflow) {
  const int2 conv = hostLongLongToInt63(val);
  *primary  = conv.x;
  *overflow = conv.y;
}

//-------------------------------------------------------------------------------------------------
void hostLongLongToInt63(const llint* val, int* primary, int* overflow, const size_t n_values) {
  for (size_t i = 0; i < n_values; i++) {
    const int2 conv = hostLongLongToInt63(val[i]);
    primary[i]  = conv.x;
    overflow[i] = conv.y;
  }
}

//-------------------------------------------------------------------------------------------------
void hostLongLongToInt63(const llint* val_x, const llint* val_y, const llint* val_z,
                         int* primary_x, int* overflow_x, int* primary_y, int* overflow_y,
                         int* primary_z, int* overflow_z, const size_t n_values) {
  hostLongLongToInt63(val_x, primary_x, overflow_x, n_values);
  hostLongLongToInt63(val_y, primary_y, overflow_y, n_values);
  hostLongLongToInt63(val_z, primary_z, overflow_z, n_values);
}

//-------------------------------------------------------------------------------------------------
void hostLongLongToInt63(const std::vector<llint> &val, std::vector<int> *primary,
                         std::vector<int> *overflow) {
  checkFPVectorLength(primary->size(), overflow->size(), "hostLongLongToInt63");
  checkFPVectorLength(primary->size(), val.size(), "hostLongLongInt63", "The original data must "
                      "have the same length as the integral result vectors.");
  hostLongLongToInt63(val.data(), primary->data(), overflow->data(), val.size());
}

//-------------------------------------------------------------------------------------------------
void hostLongLongToInt63(const std::vector<llint> &val_x, const std::vector<llint> &val_y,
                         const std::vector<llint> &val_z, std::vector<int> *primary_x,
                         std::vector<int> *overflow_x, std::vector<int> *primary_y,
                         std::vector<int> *overflow_y, std::vector<int> *primary_z,
                         std::vector<int> *overflow_z) {
  hostLongLongToInt63(val_x, primary_x, overflow_x);
  hostLongLongToInt63(val_y, primary_y, overflow_y);
  hostLongLongToInt63(val_z, primary_z, overflow_z);
}

//-------------------------------------------------------------------------------------------------
void hostLongLongToInt63(const Hybrid<llint> &val, Hybrid<int> *primary, Hybrid<int> *overflow) {
  checkFPVectorLength(primary->size(), overflow->size(), "hostLongLongToInt63");
  checkFPVectorLength(primary->size(), val.size(), "hostLongLongInt63", "The original data must "
                      "have the same length as the integral result vectors.");
  hostLongLongToInt63(val.data(), primary->data(), overflow->data(), val.size());
}

//-------------------------------------------------------------------------------------------------
void hostLongLongToInt63(const Hybrid<llint> &val_x, const Hybrid<llint> &val_y,
                         const Hybrid<llint> &val_z, Hybrid<int> *primary_x,
                         Hybrid<int> *overflow_x, Hybrid<int> *primary_y, Hybrid<int> *overflow_y,
                         Hybrid<int> *primary_z, Hybrid<int> *overflow_z) {
  hostLongLongToInt63(val_x, primary_x, overflow_x);
  hostLongLongToInt63(val_y, primary_y, overflow_y);
  hostLongLongToInt63(val_z, primary_z, overflow_z);
}

//-------------------------------------------------------------------------------------------------
int2 hostDoubleToInt63(const double dval) {
  int2 result;
  if (fabs(dval) >= max_int_accumulation) {
    const int spillover = dval / max_int_accumulation;
    result.x = dval - (static_cast<double>(spillover) * max_int_accumulation);
    result.y = spillover;
  }
  else {
    result.x = dval;
    result.y = 0;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void hostDoubleToInt63(const double dval, int *primary, int *overflow) {
  const int2 result = hostDoubleToInt63(dval);
  *primary  = result.x;
  *overflow = result.y;
}

//-------------------------------------------------------------------------------------------------
int95_t hostDoubleToInt95(const double dval) {
  int95_t result;
  if (fabs(dval) >= max_llint_accumulation) {
    const int spillover = dval / max_llint_accumulation;
    result.x = dval - (static_cast<double>(spillover) * max_llint_accumulation);
    result.y = spillover;
  }
  else {
    result.x = dval;
    result.y = 0;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void hostDoubleToInt95(const double dval, llint *primary, int *overflow) {
  const int95_t result = hostDoubleToInt95(dval);
  *primary  = result.x;
  *overflow = result.y;
}

//-------------------------------------------------------------------------------------------------
void hostDoubleToInt95(const double* dval, llint* primary, int* overflow, const size_t n_values,
                       const double scale) {
  for (size_t i = 0LLU; i < n_values; i++) {
    hostDoubleToInt95(dval[i] * scale, &primary[i], &overflow[i]);
  }
}

//-------------------------------------------------------------------------------------------------
void hostDoubleToInt95(const std::vector<double> &dval, std::vector<llint> *primary,
                       std::vector<int> *overflow, const double scale) {
  checkFPVectorLength(primary->size(), overflow->size(), "hostDoubleToInt95");
  checkFPVectorLength(primary->size(), dval.size(), "hostDoubleToInt95", "The original data must "
                      "have the same length as the integral result vectors.");
  hostDoubleToInt95(dval.data(), primary->data(), overflow->data(), dval.size(), scale);
}

//-------------------------------------------------------------------------------------------------
void hostDoubleToInt95(const Hybrid<double> &dval, Hybrid<llint> *primary, Hybrid<int> *overflow,
                   const double scale) {
  checkFPVectorLength(primary->size(), overflow->size(), "hostDoubleToInt95");
  checkFPVectorLength(primary->size(), dval.size(), "hostDoubleToInt95", "The original data must "
                      "have the same length as the integral result vectors.");
  hostDoubleToInt95(dval.data(), primary->data(), overflow->data(), dval.size(), scale);
}

//-------------------------------------------------------------------------------------------------
void hostDoubleToInt95(const double* dval_x, const double* dval_y, const double* dval_z,
                       llint* primary_x, int* overflow_x, llint* primary_y, int* overflow_y,
                       llint* primary_z, int* overflow_z, const size_t n_values,
                       const double scale) {
  for (size_t i = 0LLU; i < n_values; i++) {
    hostDoubleToInt95(dval_x[i] * scale, &primary_x[i], &overflow_x[i]);
    hostDoubleToInt95(dval_y[i] * scale, &primary_y[i], &overflow_y[i]);
    hostDoubleToInt95(dval_z[i] * scale, &primary_z[i], &overflow_z[i]);
  }
}

//-------------------------------------------------------------------------------------------------
void hostDoubleToInt95(const std::vector<double> &dval_x, const std::vector<double> &dval_y,
                       const std::vector<double> &dval_z, std::vector<llint> *primary_x,
                       std::vector<int> *overflow_x, std::vector<llint> *primary_y,
                       std::vector<int> *overflow_y, std::vector<llint> *primary_z,
                       std::vector<int> *overflow_z, const double scale) {
  hostDoubleToInt95(dval_x, primary_x, overflow_x, scale);
  hostDoubleToInt95(dval_y, primary_y, overflow_y, scale);
  hostDoubleToInt95(dval_z, primary_z, overflow_z, scale);
}

//-------------------------------------------------------------------------------------------------
void hostDoubleToInt95(const Hybrid<double> &dval_x, const Hybrid<double> &dval_y,
                       const Hybrid<double> &dval_z, Hybrid<llint> *primary_x,
                       Hybrid<int> *overflow_x, Hybrid<llint> *primary_y, Hybrid<int> *overflow_y,
                       Hybrid<llint> *primary_z, Hybrid<int> *overflow_z, const double scale) {
  hostDoubleToInt95(dval_x, primary_x, overflow_x, scale);
  hostDoubleToInt95(dval_y, primary_y, overflow_y, scale);
  hostDoubleToInt95(dval_z, primary_z, overflow_z, scale);
}

//-------------------------------------------------------------------------------------------------
llint hostInt63ToLongLong(const int primary, const int overflow) {
  llint result = overflow;
  result *= max_int_accumulation_ll;
  result += primary;
  return result;
}

//-------------------------------------------------------------------------------------------------
void hostInt63ToLongLong(llint* result, const int* primary, const int* overflow,
                         const size_t n_values) {
  for (size_t i = 0; i < n_values; i++) {
    result[i] = hostInt63ToLongLong(primary[i], overflow[i]);
  }
}

//-------------------------------------------------------------------------------------------------
void hostInt63ToLongLong(std::vector<llint> *result, const std::vector<int> &primary,
                         const std::vector<int> &overflow) {
  checkFPVectorLength(primary.size(), overflow.size(), "hostInt63ToLongLong");
  checkFPVectorLength(primary.size(), result->size(), "hostInt63ToLongLong", "The result must "
                      "have the same length as the original, integral data.");
  hostInt63ToLongLong(result->data(), primary.data(), overflow.data(), result->size());
}

//-------------------------------------------------------------------------------------------------
void hostInt63ToLongLong(Hybrid<llint> *result, const Hybrid<int> &primary,
                         const Hybrid<int> &overflow) {
  checkFPVectorLength(primary.size(), overflow.size(), "hostInt63ToLongLong");
  checkFPVectorLength(primary.size(), result->size(), "hostInt63ToLongLong", "The result must "
                      "have the same length as the original, integral data.");
  hostInt63ToLongLong(result->data(), primary.data(), overflow.data(), result->size());
}

//-------------------------------------------------------------------------------------------------
void hostInt63ToLongLong(llint* result_x, llint* result_y, llint* result_z, const int* primary_x,
                         const int* overflow_x, const int* primary_y, const int* overflow_y,
                         const int* primary_z, const int* overflow_z, const size_t n_values) {
  hostInt63ToLongLong(result_x, primary_x, overflow_x, n_values);
  hostInt63ToLongLong(result_y, primary_y, overflow_y, n_values);
  hostInt63ToLongLong(result_z, primary_z, overflow_z, n_values);
}

//-------------------------------------------------------------------------------------------------
void hostInt63ToLongLong(std::vector<llint> *result_x, std::vector<llint> *result_y,
                         std::vector<llint> *result_z, const std::vector<int> &primary_x,
                         const std::vector<int> &overflow_x, const std::vector<int> &primary_y,
                         const std::vector<int> &overflow_y, const std::vector<int> &primary_z,
                         const std::vector<int> &overflow_z) {
  hostInt63ToLongLong(result_x, primary_x, overflow_x);
  hostInt63ToLongLong(result_y, primary_y, overflow_y);
  hostInt63ToLongLong(result_z, primary_z, overflow_z);
}

//-------------------------------------------------------------------------------------------------
void hostInt63ToLongLong(Hybrid<llint> *result_x, Hybrid<llint> *result_y, Hybrid<llint> *result_z,
                         const Hybrid<int> &primary_x, const Hybrid<int> &overflow_x,
                         const Hybrid<int> &primary_y, const Hybrid<int> &overflow_y,
                         const Hybrid<int> &primary_z, const Hybrid<int> &overflow_z) {
  hostInt63ToLongLong(result_x, primary_x, overflow_x);
  hostInt63ToLongLong(result_y, primary_y, overflow_y);
  hostInt63ToLongLong(result_z, primary_z, overflow_z);
}

//-------------------------------------------------------------------------------------------------
double hostInt63ToDouble(const int primary, const int overflow) {
  return (static_cast<double>(overflow) * max_int_accumulation) + static_cast<double>(primary);
}

//-------------------------------------------------------------------------------------------------
void hostInt63ToDouble(double* result, const int* primary, const int* overflow,
                   const size_t n_values, const double descale) {
  for (size_t i = 0LLU; i < n_values; i++) {
    result[i] = hostInt63ToDouble(primary[i], overflow[i]) * descale;
  }
}

//-------------------------------------------------------------------------------------------------
void hostInt63ToDouble(std::vector<double> *result, const std::vector<int> &primary,
                   const std::vector<int> &overflow, const double descale) {
  checkFPVectorLength(primary.size(), overflow.size(), "hostInt63ToDouble");
  checkFPVectorLength(primary.size(), result->size(), "hostInt63ToDouble", "The result must "
                      "have the same length as the original, integral data.");
  hostInt63ToDouble(result->data(), primary.data(), overflow.data(), primary.size(), descale);
}

//-------------------------------------------------------------------------------------------------
void hostInt63ToDouble(Hybrid<double> *result, const Hybrid<int> &primary,
                       const Hybrid<int> &overflow, const double descale) {
  checkFPVectorLength(primary.size(), overflow.size(), "hostInt63ToDouble");
  checkFPVectorLength(primary.size(), result->size(), "hostInt63ToDouble", "The result must "
                      "have the same length as the original, integral data.");
  hostInt63ToDouble(result->data(), primary.data(), overflow.data(), primary.size(), descale);
}

//-------------------------------------------------------------------------------------------------
float hostInt63ToFloat(const int primary, const int overflow) {
  return (static_cast<float>(overflow) * max_int_accumulation_f) + static_cast<float>(primary);
}

//-------------------------------------------------------------------------------------------------
void hostInt63ToFloat(float* result, const int* primary, const int* overflow,
                      const size_t n_values, const float descale) {
  for (size_t i = 0LLU; i < n_values; i++) {
    result[i] = hostInt63ToFloat(primary[i], overflow[i]) * descale;
  }
}

//-------------------------------------------------------------------------------------------------
void hostInt63ToFloat(std::vector<float> *result, const std::vector<int> &primary,
                  const std::vector<int> &overflow, const float descale) {
  checkFPVectorLength(primary.size(), overflow.size(), "hostInt63ToFloat");
  checkFPVectorLength(primary.size(), result->size(), "hostInt63ToFloat", "The result must "
                      "have the same length as the original, integral data.");
  hostInt63ToFloat(result->data(), primary.data(), overflow.data(), primary.size(), descale);
}

//-------------------------------------------------------------------------------------------------
void hostInt63ToFloat(Hybrid<float> *result, const Hybrid<int> &primary,
                      const Hybrid<int> &overflow, const float descale) {
  checkFPVectorLength(primary.size(), overflow.size(), "hostInt63ToFloat");
  checkFPVectorLength(primary.size(), result->size(), "hostInt63ToFloat", "The result must "
                      "have the same length as the original, integral data.");
  hostInt63ToFloat(result->data(), primary.data(), overflow.data(), primary.size(), descale);
}

//-------------------------------------------------------------------------------------------------
void hostInt63ToDouble(double* result_x, double* result_y, double* result_z, const int* primary_x,
                       const int* overflow_x, const int* primary_y, const int* overflow_y,
                       const int* primary_z, const int* overflow_z, const size_t n_values,
                       const double descale) {
  for (size_t i = 0LLU; i < n_values; i++) {
    result_x[i] = hostInt63ToDouble(primary_x[i], overflow_x[i]) * descale;
    result_y[i] = hostInt63ToDouble(primary_y[i], overflow_y[i]) * descale;
    result_z[i] = hostInt63ToDouble(primary_z[i], overflow_z[i]) * descale;
  }
}

//-------------------------------------------------------------------------------------------------
void hostInt63ToDouble(std::vector<double> *result_x, std::vector<double> *result_y,
                       std::vector<double> *result_z, const std::vector<int> &primary_x,
                       const std::vector<int> &overflow_x, const std::vector<int> &primary_y,
                       const std::vector<int> &overflow_y, const std::vector<int> &primary_z,
                       const std::vector<int> &overflow_z, const size_t n_values,
                       const double descale) {
  hostInt63ToDouble(result_x, primary_x, overflow_x, descale);
  hostInt63ToDouble(result_y, primary_y, overflow_y, descale);
  hostInt63ToDouble(result_z, primary_z, overflow_z, descale);
}

//-------------------------------------------------------------------------------------------------
void hostInt63ToDouble(Hybrid<double> *result_x, Hybrid<double> *result_y,
                       Hybrid<double> *result_z, const Hybrid<int> &primary_x,
                       const Hybrid<int> &overflow_x, const Hybrid<int> &primary_y,
                       const Hybrid<int> &overflow_y, const Hybrid<int> &primary_z,
                       const Hybrid<int> &overflow_z, const size_t n_values,
                       const double descale) {
  hostInt63ToDouble(result_x, primary_x, overflow_x, descale);
  hostInt63ToDouble(result_y, primary_y, overflow_y, descale);
  hostInt63ToDouble(result_z, primary_z, overflow_z, descale);
}

//-------------------------------------------------------------------------------------------------
void hostInt63ToFloat(float* result_x, float* result_y, float* result_z, const int* primary_x,
                      const int* overflow_x, const int* primary_y, const int* overflow_y,
                      const int* primary_z, const int* overflow_z, const size_t n_values,
                      const float descale) {
  for (size_t i = 0LLU; i < n_values; i++) {
    result_x[i] = hostInt63ToFloat(primary_x[i], overflow_x[i]) * descale;
    result_y[i] = hostInt63ToFloat(primary_y[i], overflow_y[i]) * descale;
    result_z[i] = hostInt63ToFloat(primary_z[i], overflow_z[i]) * descale;
  }
}

//-------------------------------------------------------------------------------------------------
void hostInt63ToFloat(std::vector<float> *result_x, std::vector<float> *result_y,
                      std::vector<float> *result_z, const std::vector<int> &primary_x,
                      const std::vector<int> &overflow_x, const std::vector<int> &primary_y,
                      const std::vector<int> &overflow_y, const std::vector<int> &primary_z,
                      const std::vector<int> &overflow_z, const float descale) {
  hostInt63ToFloat(result_x, primary_x, overflow_x, descale);
  hostInt63ToFloat(result_y, primary_y, overflow_y, descale);
  hostInt63ToFloat(result_z, primary_z, overflow_z, descale);
}

//-------------------------------------------------------------------------------------------------
void hostInt63ToFloat(Hybrid<float> *result_x, Hybrid<float> *result_y, Hybrid<float> *result_z,
                      const Hybrid<int> &primary_x, const Hybrid<int> &overflow_x,
                      const Hybrid<int> &primary_y, const Hybrid<int> &overflow_y,
                      const Hybrid<int> &primary_z, const Hybrid<int> &overflow_z,
                      const float descale) {
  hostInt63ToFloat(result_x, primary_x, overflow_x, descale);
  hostInt63ToFloat(result_y, primary_y, overflow_y, descale);
  hostInt63ToFloat(result_z, primary_z, overflow_z, descale);
}

//-------------------------------------------------------------------------------------------------
double hostInt63ToDouble(const int2 ival) {
  return (static_cast<double>(ival.y) * max_int_accumulation) + static_cast<double>(ival.x);
}

//-------------------------------------------------------------------------------------------------
float hostInt63ToFloat(const int2 ival) {
  return (static_cast<float>(ival.y) * max_int_accumulation_f) + static_cast<float>(ival.x);
}

//-------------------------------------------------------------------------------------------------
double hostInt95ToDouble(const int95_t ival) {
  return (static_cast<double>(ival.y) * max_llint_accumulation) + static_cast<double>(ival.x);
}

//-------------------------------------------------------------------------------------------------
double hostInt95ToDouble(const llint primary, const int overflow) {
  return (static_cast<double>(overflow) * max_llint_accumulation) + static_cast<double>(primary);
}

//-------------------------------------------------------------------------------------------------
void hostInt95ToDouble(double* result, const llint* primary, const int* overflow,
                   const size_t n_values, const double descale) {
  for (size_t i = 0LLU; i < n_values; i++) {
    result[i] = hostInt95ToDouble(primary[i], overflow[i]) * descale;
  }
}

//-------------------------------------------------------------------------------------------------
void hostInt95ToDouble(std::vector<double> *result, const std::vector<llint> &primary,
                   const std::vector<int> &overflow, const double descale) {
  checkFPVectorLength(primary.size(), overflow.size(), "hostInt95ToDouble");
  checkFPVectorLength(primary.size(), result->size(), "hostInt95ToDouble", "The result must "
                      "have the same length as the original, integral data.");
  hostInt95ToDouble(result->data(), primary.data(), overflow.data(), primary.size(), descale);
}

//-------------------------------------------------------------------------------------------------
void hostInt95ToDouble(Hybrid<double> *result, const Hybrid<llint> &primary,
                   const Hybrid<int> &overflow, const double descale) {
  checkFPVectorLength(primary.size(), overflow.size(), "hostInt95ToDouble");
  checkFPVectorLength(primary.size(), result->size(), "hostInt95ToDouble", "The result must "
                      "have the same length as the original, integral data.");
  hostInt95ToDouble(result->data(), primary.data(), overflow.data(), primary.size(), descale);
}

//-------------------------------------------------------------------------------------------------
void hostInt95ToDouble(double* result_x, double* result_y, double* result_z,
                       const llint* primary_x, const int* overflow_x, const llint* primary_y,
                       const int* overflow_y, const llint* primary_z, const int* overflow_z,
                       const size_t n_values, const double descale) {
  for (size_t i = 0LLU; i < n_values; i++) {
    result_x[i] = hostInt95ToDouble(primary_x[i], overflow_x[i]) * descale;
    result_y[i] = hostInt95ToDouble(primary_y[i], overflow_y[i]) * descale;
    result_z[i] = hostInt95ToDouble(primary_z[i], overflow_z[i]) * descale;
  }
}

//-------------------------------------------------------------------------------------------------
void hostInt95ToDouble(std::vector<double> *result_x, std::vector<double> *result_y,
                   std::vector<double> *result_z, const std::vector<llint> &primary_x,
                   const std::vector<int> &overflow_x, const std::vector<llint> &primary_y,
                   const std::vector<int> &overflow_y, const std::vector<llint> &primary_z,
                   const std::vector<int> &overflow_z, const double descale) {
  hostInt95ToDouble(result_x, primary_x, overflow_x, descale);
  hostInt95ToDouble(result_y, primary_y, overflow_y, descale);
  hostInt95ToDouble(result_z, primary_z, overflow_z, descale);
}

//-------------------------------------------------------------------------------------------------
void hostInt95ToDouble(Hybrid<double> *result_x, Hybrid<double> *result_y,
                       Hybrid<double> *result_z, const Hybrid<llint> &primary_x,
                       const Hybrid<int> &overflow_x, const Hybrid<llint> &primary_y,
                       const Hybrid<int> &overflow_y, const Hybrid<llint> &primary_z,
                       const Hybrid<int> &overflow_z, const double descale) {
  hostInt95ToDouble(result_x, primary_x, overflow_x, descale);
  hostInt95ToDouble(result_y, primary_y, overflow_y, descale);
  hostInt95ToDouble(result_z, primary_z, overflow_z, descale);
}

//-------------------------------------------------------------------------------------------------
void hostSplitAccumulation(const double dval, llint *primary, int *overflow) {
  llint ival;
  if (fabs(dval) >= max_llint_accumulation) {
    const int spillover = dval / max_llint_accumulation;
    ival = dval - (static_cast<double>(spillover) * max_llint_accumulation);
    *overflow += spillover;
  }
  else {
    ival = dval;
  }
  const llint prim_old = *primary;
  *primary += ival;
  const llint prim_old_plus_ival = prim_old + ival;
  if ((prim_old ^ prim_old_plus_ival) < 0LL && (prim_old ^ ival) >= 0LL) {
    *overflow += (1 - (2 * (ival < 0LL))) * 2;
  }
}

//-------------------------------------------------------------------------------------------------
void hostSplitAccumulation(const float fval, int *primary, int *overflow) {
  int ival;
  if (fabsf(fval) >= max_int_accumulation_f) {
    const int spillover = fval / max_int_accumulation_f;
    ival = fval - (static_cast<float>(spillover) * max_int_accumulation_f);
    *overflow += spillover;
  }
  else {
    ival = fval;
  }
  const int prim_old = *primary;
  *primary += ival;
  const int prim_old_plus_ival = prim_old + ival;
  if ((prim_old ^ prim_old_plus_ival) < 0 && (prim_old ^ ival) >= 0) {
    *overflow += (1 - (2 * (ival < 0))) * 2;
  }
}

//-------------------------------------------------------------------------------------------------
int95_t hostSplitFPSum(const int95_t a, const int95_t b) {
  int95_t result = { a.x + b.x, a.y + b.y };
  result.y += (1 - (2 * (b.x < 0LL))) * ((a.x ^ result.x) < 0 && (a.x ^ b.x) >= 0LL) * 2;
  return result;
}

//-------------------------------------------------------------------------------------------------
int2 hostSplitFPSum(const int2 a, const int2 b) {
  int2 result = { a.x + b.x, a.y + b.y };
  result.y += (1 - (2 * (b.x < 0))) * ((a.x ^ result.x) < 0 && (a.x ^ b.x) >= 0) * 2;
  return result;
}

//-------------------------------------------------------------------------------------------------
int95_t hostSplitFPSum(const int95_t a, const double breal) {
  const int95_t b = hostDoubleToInt95(breal);
  int95_t result = { a.x + b.x, a.y + b.y };
  result.y += (1 - (2 * (b.x < 0LL))) * ((a.x ^ result.x) < 0 && (a.x ^ b.x) >= 0LL) * 2;
  return result;
}

//-------------------------------------------------------------------------------------------------
int2 hostSplitFPSum(const int2 a, const float breal) {
  const int2 b = hostFloatToInt63(breal);
  int2 result = { a.x + b.x, a.y + b.y };
  result.y += (1 - (2 * (b.x < 0))) * ((a.x ^ result.x) < 0 && (a.x ^ b.x) >= 0) * 2;
  return result;
}

//-------------------------------------------------------------------------------------------------
int95_t hostSplitFPSum(const int95_t a, const llint b_x, const int b_y) {
  int95_t result = { a.x + b_x, a.y + b_y };
  result.y += (1 - (2 * (b_x < 0LL))) * ((a.x ^ result.x) < 0 && (a.x ^ b_x) >= 0LL) * 2;
  return result;
}

//-------------------------------------------------------------------------------------------------
int2 hostSplitFPSum(const int2 a, const int b_x, const int b_y) {
  int2 result = { a.x + b_x, a.y + b_y };
  result.y += (1 - (2 * (b_x < 0))) * ((a.x ^ result.x) < 0 && (a.x ^ b_x) >= 0) * 2;
  return result;
}

//-------------------------------------------------------------------------------------------------
int95_t hostInt95Sum(const llint a_x, const int a_y, const llint b_x, const int b_y) {
  int95_t result = { a_x + b_x, a_y + b_y };
  result.y += (1 - (2 * (b_x < 0LL))) * ((a_x ^ result.x) < 0 && (a_x ^ b_x) >= 0LL) * 2;
  return result;
}

//-------------------------------------------------------------------------------------------------
int2 hostInt63Sum(const int a_x, const int a_y, const int b_x, const int b_y) {
  int2 result = { a_x + b_x, a_y + b_y };
  result.y += (1 - (2 * (b_x < 0))) * ((a_x ^ result.x) < 0 && (a_x ^ b_x) >= 0) * 2;
  return result;
}

//-------------------------------------------------------------------------------------------------
int95_t hostInt95Sum(const llint a_x, const int a_y, const double breal) {
  const int95_t b = hostDoubleToInt95(breal);
  int95_t result = { a_x + b.x, a_y + b.y };
  result.y += (1 - (2 * (b.x < 0LL))) * ((a_x ^ result.x) < 0 && (a_x ^ b.x) >= 0LL) * 2;
  return result;
}

//-------------------------------------------------------------------------------------------------
int2 hostInt63Sum(const int a_x, const int a_y, const float breal) {
  const int2 b = hostFloatToInt63(breal);
  int2 result = { a_x + b.x, a_y + b.y };
  result.y += (1 - (2 * (b.x < 0))) * ((a_x ^ result.x) < 0 && (a_x ^ b.x) >= 0) * 2;
  return result;
}

//-------------------------------------------------------------------------------------------------
int95_t hostSplitFPSubtract(const int95_t a, const int95_t b) {
  const int95_t neg_b = { -b.x, -b.y + (2 * (b.x == LLONG_MIN)) };
  return hostSplitFPSum(a, neg_b);
}

//-------------------------------------------------------------------------------------------------
int2 hostSplitFPSubtract(const int2 a, const int2 b) {
  const int2 neg_b = { -b.x, -b.y + (2 * (b.x == INT_MIN)) };
  return hostSplitFPSum(a, neg_b);
}

//-------------------------------------------------------------------------------------------------
int95_t hostSplitFPSubtract(const int95_t a, const llint b_x, const int b_y) {
  const int95_t neg_b = { -b_x, -b_y + (2 * (b_x == LLONG_MIN)) };
  return hostSplitFPSum(a, neg_b);
}

//-------------------------------------------------------------------------------------------------
int2 hostSplitFPSubtract(const int2 a, const int b_x, const int b_y) {
  const int2 neg_b = { -b_x, -b_y + (2 * (b_x == INT_MIN)) };
  return hostSplitFPSum(a, neg_b);
}

//-------------------------------------------------------------------------------------------------
int95_t hostInt95Subtract(const llint a_x, const int a_y, const llint b_x, const int b_y) {
  return hostInt95Sum(a_x, a_y, -b_x, -b_y + (2 * (b_x == LLONG_MIN)));
}

//-------------------------------------------------------------------------------------------------
int2 hostInt63Subtract(const int a_x, const int a_y, const int b_x, const int b_y) {
  return hostInt63Sum(a_x, a_y, -b_x, -b_y + (2 * (b_x == INT_MIN)));
}

//-------------------------------------------------------------------------------------------------
int95_t hostSplitFPMult(const int95_t a, const int b) {

  // Break the split fixed-precision number into two halves of 47 high bits and 48 low bits.
  const llint ll_b = b;
  llint low = (a.x & 0xffffffffLL) * ll_b;
  llint mid = ((static_cast<ullint>(a.x) & 0xffffffff00000000LLU) >> 32);
  if (a.x < 0LL) {
    mid |= 0xffffffff00000000LL;
  }

  // If a.x is the most negative possible long long integer and b is a very large number, the
  // format will be broken.  Such an event is unavoidable and will not be trapped, i.e. if mid
  // times b breaks the long long integer format, the int95_t format would be broken.  Similarly,
  // if a.y times b breaks the int32_t format, that will break int95_t and likewise needs not be
  // trapped.
  mid *= ll_b;

  // The mid value measures how many portions of 2^32 that the product now contains.  The overflow
  // of the result will measure how many portions of 2^3 the product contains, to which mid will
  // contribute its value divided by 2^31 (that is, max_int_accumulation as defined in the library
  // header split_fixed_precision.h).  The low value, even if multiplied by the largest possible
  // int32_t value, will not contribute to the overflow bits in the result. 
  const int mid_ycontrib = mid / max_int_accumulation_ll;
  const llint mid_xcontrib = (mid - (max_int_accumulation_ll * static_cast<llint>(mid_ycontrib))) *
                             (0x0000000100000000LL);
  return hostInt95Sum(mid_xcontrib, mid_ycontrib, low, a.y * b);
}

//-------------------------------------------------------------------------------------------------
int95_t hostInt95Mult(const llint a_x, const int a_y, const int b) {
  const int95_t a = { a_x, a_y };
  return hostSplitFPMult(a, b);
}

//-------------------------------------------------------------------------------------------------
int2 hostSplitFPMult(const int2 a, const int b) {
  const llint intermediate = hostInt63ToLongLong(a.x, a.y);
  return hostLongLongToInt63(intermediate * b);
}

//-------------------------------------------------------------------------------------------------
int2 hostInt63Mult(const int a_x, const int a_y, const int b) {
  const int2 a = { a_x, a_y };
  return hostSplitFPMult(a, b);
}

//-------------------------------------------------------------------------------------------------
int2 hostChangeFPBits(const int2 fp, const int native_bits, const int output_bits) {
  if (native_bits == output_bits) {
    return fp;
  }

  // Compute the conversion factor, relying on the fact that a double (64-bit floating point
  // number) will align with a long long unsigned int.
  PolyNumeric conv_factor;
  conv_factor.ulli = 1023 + output_bits - native_bits;
  conv_factor.ulli <<= 52;
  const double xcomp = static_cast<double>(fp.x) * conv_factor.d;
  const double ycomp = static_cast<double>(fp.y) * max_int_accumulation * conv_factor.d;
  const int2 xnew = hostDoubleToInt63(xcomp);
  const int2 ynew = hostDoubleToInt63(ycomp);
  return hostSplitFPSum(xnew, ynew);  
}
  
//-------------------------------------------------------------------------------------------------
int95_t hostChangeFPBits(const int95_t fp, const int native_bits, const int output_bits) {
  if (native_bits == output_bits) {
    return fp;
  }

  // Break the 64-bit floating point number into two 32-bit components, each of which can be
  // handled exactly by a 64-bit floating point number.
  PolyNumeric conv_factor;
  conv_factor.ulli = 1023 + output_bits - native_bits;
  conv_factor.ulli <<= 52;
  llint ilow_xcomp = (static_cast<ullint>(fp.x) & 0xffffffffLLU);
  llint ihigh_xcomp = fp.x - ilow_xcomp;
  const double xcomp_low  = static_cast<double>(ilow_xcomp) * conv_factor.d;
  const double xcomp_high = static_cast<double>(ihigh_xcomp) * conv_factor.d;
  const double ycomp = static_cast<double>(fp.y) * max_llint_accumulation * conv_factor.d;
  const int95_t xnew_low  = hostDoubleToInt95(xcomp_low);
  const int95_t xnew_high = hostDoubleToInt95(xcomp_high);
  const int95_t xnew      = hostSplitFPSum(xnew_low, xnew_high);
  const int95_t ynew      = hostDoubleToInt95(ycomp);
  return hostSplitFPSum(xnew, ynew);
}

//-------------------------------------------------------------------------------------------------
void hostChangeFPBits(std::vector<int95_t> *fp, const int native_bits, const int output_bits) {
  int95_t* fp_ptr = fp->data();
  const size_t length = fp->size();
  for (size_t i = 0; i < length; i++) {
    fp_ptr[i] = hostChangeFPBits(fp_ptr[i], native_bits, output_bits);
  }
}

//-------------------------------------------------------------------------------------------------
void hostChangeFPBits(int* fp, int* fp_ovrf, const size_t length, const int native_bits,
                      const int output_bits) {
  for (size_t i = 0; i < length; i++) {
    const int2 orig = { fp[i], fp_ovrf[i] };
    const int2 conv = hostChangeFPBits(orig, native_bits, output_bits);
    fp[i] = conv.x;
    fp_ovrf[i] = conv.y;
  }
}

//-------------------------------------------------------------------------------------------------
void hostChangeFPBits(llint* fp, int* fp_ovrf, const size_t length, const int native_bits,
                      const int output_bits) {
  for (size_t i = 0; i < length; i++) {
    const int95_t orig = { fp[i], fp_ovrf[i] };
    const int95_t conv = hostChangeFPBits(orig, native_bits, output_bits);
    fp[i] = conv.x;
    fp_ovrf[i] = conv.y;
  }
}

//-------------------------------------------------------------------------------------------------
void hostChangeFPBits(std::vector<int> *fp, std::vector<int> *fp_ovrf, const int native_bits,
                      const int output_bits) {
  checkFPVectorLength(fp->size(), fp_ovrf->size(), "hostChangeFPBits");
  hostChangeFPBits(fp->data(), fp_ovrf->data(), fp->size(), native_bits, output_bits);
}

//-------------------------------------------------------------------------------------------------
void hostChangeFPBits(std::vector<llint> *fp, std::vector<int> *fp_ovrf, const int native_bits,
                      const int output_bits) {
  checkFPVectorLength(fp->size(), fp_ovrf->size(), "hostChangeFPBits");
  hostChangeFPBits(fp->data(), fp_ovrf->data(), fp->size(), native_bits, output_bits);
}

//-------------------------------------------------------------------------------------------------
void hostChangeFPBits(Hybrid<int> *fp, Hybrid<int> *fp_ovrf, const int native_bits,
                      const int output_bits) {
  checkFPVectorLength(fp->size(), fp_ovrf->size(), "hostChangeFPBits");
  hostChangeFPBits(fp->data(), fp_ovrf->data(), fp->size(), native_bits, output_bits);
}

//-------------------------------------------------------------------------------------------------
void hostChangeFPBits(Hybrid<llint> *fp, Hybrid<int> *fp_ovrf, const int native_bits,
                      const int output_bits) {
  checkFPVectorLength(fp->size(), fp_ovrf->size(), "hostChangeFPBits");
  hostChangeFPBits(fp->data(), fp_ovrf->data(), fp->size(), native_bits, output_bits);
}

//-------------------------------------------------------------------------------------------------
void fixedPrecisionGrid(std::vector<int95_t> *coordinates, const int95_t origin,
                        const int95_t increment) {
  const size_t grid_size = coordinates->size();
  int95_t marker = origin;
  int95_t* coord_ptr = coordinates->data();
  for (size_t i = 0LLU; i < grid_size; i++) {
    coord_ptr[i] = marker;
    marker = hostSplitFPSum(marker, increment);
  }
}

//-------------------------------------------------------------------------------------------------
void fixedPrecisionGrid(std::vector<llint> *primary, std::vector<int> *overflow,
                        const int95_t origin, const int95_t increment) {
  checkFPVectorLength(primary->size(), overflow->size(), "hostChangeFPBits", "Both primary and "
                      "overflow vectors must be pre-allocated to the same size.  Lengths "
                      "are " + std::to_string(primary->size()) + " and " +
                      std::to_string(overflow->size()) + ".");
  fixedPrecisionGrid(primary->data(), overflow->data(), origin, increment, primary->size());
}

//-------------------------------------------------------------------------------------------------
void fixedPrecisionGrid(Hybrid<llint> *primary, Hybrid<int> *overflow,
                        const int95_t origin, const int95_t increment) {
  checkFPVectorLength(primary->size(), overflow->size(), "hostChangeFPBits", "Both primary and "
                      "overflow vectors must be pre-allocated to the same size.  Lengths "
                      "are " + std::to_string(primary->size()) + " and " +
                      std::to_string(overflow->size()) + ".");
  fixedPrecisionGrid(primary->data(), overflow->data(), origin, increment, primary->size());
}

//-------------------------------------------------------------------------------------------------
void fixedPrecisionGrid(llint *primary, int *overflow, const int95_t origin,
                        const int95_t increment, const size_t grid_size) {
  int95_t marker = origin;
  for (size_t i = 0LLU; i < grid_size; i++) {
    primary[i] = marker.x;
    overflow[i] = marker.y;
    marker = hostSplitFPSum(marker, increment);
  }
}
  
} // namespace numerics
} // namespace stormm
