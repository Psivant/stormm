// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
template <typename T> void vectorComparisonCheck(const std::vector<T> &va,
                                                 const std::vector<T> &vb, const char* caller) {

  // Trap comparisons of non-scalar types
  if (isScalarType<T>() == false) {
    rtErr("Comparison between vectors of type " + std::string(typeid(T).name()) +
          " is not permitted.", caller);
  }

  // Trap vector size mismatch
  if (va.size() != vb.size()) {
    rtErr("Comparison requires vectors of identical sizes (" + std::to_string(va.size()) +
          " and " + std::to_string(vb.size()) + " given).", caller);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void vectorComparisonCheck(const Hybrid<T> &va, const Hybrid<T> &vb,
                                                 const char* caller) {

  // Trap comparisons of non-scalar types
  if (isScalarType<T>() == false) {
    rtErr("Comparison between vectors of type " + std::string(typeid(T).name()) +
          " is not permitted.", caller);
  }

  // Trap vector size mismatch
  if (va.size() != vb.size()) {
    rtErr("Comparison requires vectors of identical sizes (" + std::to_string(va.size()) +
          " and " + std::to_string(vb.size()) + " given).", caller);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> T maxAbsoluteDifference(const T* va, const T* vb, const size_t length) {
  if (isSignedIntegralScalarType<T>()) {
    T mdev = 0;
    for (size_t i = 0; i < length; i++) {
      T tdev = abs(va[i] - vb[i]);
      if (tdev > mdev) {
        mdev = tdev;
      }
    }
    return mdev;
  }
  else if (isUnsignedIntegralScalarType<T>()) {
    T mdev = 0;
    for (size_t i = 0; i < length; i++) {
      T tdev = (va[i] > vb[i]) ? va[i] - vb[i] : vb[i] - va[i];
      if (tdev > mdev) {
        mdev = tdev;
      }
    }
    return mdev;
  }
  else if (isFloatingPointScalarType<T>()) {
    T mdev = 0;
    for (size_t i = 0; i < length; i++) {
      T tdev = fabs(va[i] - vb[i]);
      if (tdev > mdev) {
        mdev = tdev;
      }
    }
    return mdev;
  }
  else {
    rtErr("Data of type " + std::string(typeid(T).name()) + " is not suitable for comparisons.  "
          "Update isScalarType() and its related subclassification functions.",
          "maxAbsoluteDifference");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T maxAbsoluteDifference(const std::vector<T> &va, const std::vector<T> &vb) {
  vectorComparisonCheck(va, vb, "maxAbsoluteDifference");
  return maxAbsoluteDifference(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> T maxAbsoluteDifference(const Hybrid<T> &va, const Hybrid<T> &vb) {
  vectorComparisonCheck(va, vb, "maxAbsoluteDifference");
  return maxAbsoluteDifference(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double maxRelativeDifference(const T* va, const T* vb, const size_t length) {
  double mdev = 0.0;
  double abs_mdev = 0.0;
  for (size_t i = 0; i < length; i++) {
    double dvai = static_cast<double>(va[i]);
    double dvbi = static_cast<double>(vb[i]);
    double tdev = (dvai - dvbi) / dvbi;
    double abs_tdev = fabs(tdev);
    if (abs_tdev > abs_mdev) {
      mdev = tdev;
      abs_mdev = abs_tdev;
    }
  }
  return mdev;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double maxRelativeDifference(const std::vector<T> &va,
                                                   const std::vector<T> &vb) {

  // Compute the maximum relative deviation in double precision
  vectorComparisonCheck(va, vb, "maxRelativeDifference");
  return maxRelativeDifference(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double maxRelativeDifference(const Hybrid<T> &va, const Hybrid<T> &vb) {

  // Compute the maximum relative deviation in double precision
  vectorComparisonCheck(va, vb, "maxRelativeDifference");
  return maxRelativeDifference(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double meanUnsignedError(const T* va, const T* vb, const size_t length) {

  // Accumulate the mean unsigned error in double precision
  // after computing differences in the vectors' native type
  double mue = 0.0;
  if (isSignedIntegralScalarType<T>() || isFloatingPointScalarType<T>()) {
    for (size_t i = 0; i < length; i++) {
      T tdev = va[i] - vb[i];
      mue += fabs(static_cast<double>(tdev));
    }
  }
  else if (isUnsignedIntegralScalarType<T>()) {
    for (size_t i = 0; i < length; i++) {
      T tdev = (va[i] > vb[i]) ? va[i] - vb[i] : vb[i] - va[i];
      mue += static_cast<double>(tdev);
    }
  }
  else {
    rtErr("Data of type " + std::string(typeid(T).name()) + " is not suitable for comparisons.  "
          "Update isScalarType() and its related subclassification functions.",
          "meanUnsignedError");
  }
  return mue / static_cast<double>(length);
}

//-------------------------------------------------------------------------------------------------
template <typename T> double meanUnsignedError(const std::vector<T> &va,
                                               const std::vector<T> &vb) {
  vectorComparisonCheck(va, vb, "meanUnsignedError");
  return meanUnsignedError(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double meanUnsignedError(const Hybrid<T> &va, const Hybrid<T> &vb) {
  vectorComparisonCheck(va, vb, "meanUnsignedError");
  return meanUnsignedError(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double rmsError(const T* va, const T* vb, const size_t length) {
  double rmse = 0.0;
  for (size_t i = 0; i < length; i++) {
    double dvai = static_cast<double>(va[i]);
    double dvbi = static_cast<double>(vb[i]);
    double tdev = (dvai - dvbi);
    rmse += tdev * tdev;
  }
  return sqrt(rmse / static_cast<double>(length));
}
  
//-------------------------------------------------------------------------------------------------
template <typename T> double rmsError(const std::vector<T> &va, const std::vector<T> &vb) {
  vectorComparisonCheck(va, vb, "rmsError");
  return rmsError(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double rmsError(const Hybrid<T> &va, const Hybrid<T> &vb) {
  vectorComparisonCheck(va, vb, "relativeRmsError");
  return rmsError(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double relativeRmsError(const T* va, const T* vb, const size_t length) {

  // Accumulate the root mean squared error and mean unsigned
  // value of the reference vector, all in double precision
  double rmse = 0.0;
  double mvalue = 0.0;
  for (size_t i = 0; i < length; i++) {
    double dvai = static_cast<double>(va[i]);
    double dvbi = static_cast<double>(vb[i]);
    double tdev = (dvai - dvbi);
    rmse += tdev * tdev;
    mvalue += fabs(dvbi);
  }
  rmse = sqrt(rmse / static_cast<double>(length));
  return rmse / mvalue;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double relativeRmsError(const std::vector<T> &va,
                                              const std::vector<T> &vb) {
  vectorComparisonCheck(va, vb, "relativeRmsError");
  return relativeRmsError(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double relativeRmsError(const Hybrid<T> &va, const Hybrid<T> &vb) {
  vectorComparisonCheck(va, vb, "relativeRmsError");
  return relativeRmsError(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double mean(const T* va, const size_t length) {
  const double inv_length = 1.0 / static_cast<double>(length);
  double mvalue = 0.0;
  if (isScalarType<T>()) {
    for (size_t i = 0; i < length; i++) {
      mvalue += static_cast<double>(va[i]);
    }
  }
  else {
    rtErr("Data type " + std::string(typeid(T).name()) + " is not suitable for computing a mean "
          "of scalar values.", "mean");
  }
  return mvalue * inv_length;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double mean(const std::vector<T> &va) {
  return mean(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double mean(const Hybrid<T> &va) {
  return mean(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double variance(const T* va, const size_t length,
                                      const VarianceMethod method) {
  if (isScalarType<T>() == false) {
    rtErr("Data type " + std::string(typeid(T).name()) + " is not suitable for computing the "
          "variance of scalar values.", "mean");
  }
  if (length == 0) {
    rtErr("No data was present for computation.", "variance");
  }
  double s1 = 0.0;
  double s2 = 0.0;
  double mvalue = 0.0;
  const double dlength = static_cast<double>(length);
  switch (method) {
  case VarianceMethod::VARIANCE:
  case VarianceMethod::STANDARD_DEVIATION:
  case VarianceMethod::ROOT_MEAN_SQUARED_DEVIATION:
    for (size_t i = 0; i < length; i++) {
      double tval = static_cast<double>(va[i]);
      s1 += tval;
      s2 += tval * tval;
    }
    break;
  case VarianceMethod::COEFFICIENT_OF_VARIATION:
  case VarianceMethod::NORMALIZED_RMSD:
    for (size_t i = 0; i < length; i++) {
      double tval = static_cast<double>(va[i]);
      s1 += tval;
      s2 += tval * tval;
      mvalue += tval;
    }
    mvalue /= dlength;
    break;
  }
  switch (method) {
  case VarianceMethod::VARIANCE:
    return ((dlength * s2) - (s1 * s1)) / (dlength * dlength);
  case VarianceMethod::STANDARD_DEVIATION:
    {
      if (length == 1) {
        rtErr("Standard deviation is undefined for a single number.");
      }
      const double sqrt_arg = ((dlength * s2) - (s1 * s1)) / (dlength * (dlength - 1.0));
      return (sqrt_arg < 0.0) ? 0.0 : sqrt(sqrt_arg);
    }
    break;
  case VarianceMethod::ROOT_MEAN_SQUARED_DEVIATION:
    {
      const double sqrt_arg = (dlength * s2) - (s1 * s1);
      return (sqrt_arg < 0.0) ? 0.0 : sqrt(sqrt_arg) / dlength;
    }
    break;
  case VarianceMethod::COEFFICIENT_OF_VARIATION:
    {
      const double sqrt_arg = ((dlength * s2) - (s1 * s1)) / (dlength * (dlength - 1.0));
      return (sqrt_arg < 0.0) ? 0.0 : sqrt(sqrt_arg) / fabs(mvalue);
    }
    break;
  case VarianceMethod::NORMALIZED_RMSD:
    {
      const double sqrt_arg = (dlength * s2) - (s1 * s1);
      return (sqrt_arg < 0.0) ? 0.0 : sqrt(sqrt_arg) / fabs(dlength * mvalue);
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> double variance(const std::vector<T> &va, const VarianceMethod method) {
  return variance(va.data(), va.size(), method);
}

//-------------------------------------------------------------------------------------------------
template <typename T> double variance(const Hybrid<T> &va, const VarianceMethod method) {
  return variance(va.data(), va.size(), method);
}

//-------------------------------------------------------------------------------------------------
template <typename T> double2 bivariateMean(const T* va, const size_t length) {
  if (isHpcVectorType<T>() && getHpcVectorTypeSize<T>() == 2) {
    const double inv_length = 1.0 / static_cast<double>(length);
    double2 mvalue = { 0.0, 0.0 };
    for (size_t i = 0; i < length; i++) {
      mvalue.x += static_cast<double>(va[i].x);
      mvalue.y += static_cast<double>(va[i].y);
    }
    mvalue.x *= inv_length;
    mvalue.y *= inv_length;
    return mvalue;
  }
  else {
    rtErr("Data type " + std::string(typeid(T).name()) + " is not a suitable 2-tuple.",
          "bivariateMean");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> double2 bivariateMean(const std::vector<T> &va) {
  return bivariateMean(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double2 bivariateMean(const Hybrid<T> &va) {
  return bivariateMean(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double3 trivariateMean(const T* va, const size_t length) {
  if (isHpcVectorType<T>() && getHpcVectorTypeSize<T>() == 3) {
    const double inv_length = 1.0 / static_cast<double>(length);
    double3 mvalue = { 0.0, 0.0, 0.0 };
    for (size_t i = 0; i < length; i++) {
      mvalue.x += static_cast<double>(va[i].x);
      mvalue.y += static_cast<double>(va[i].y);
      mvalue.z += static_cast<double>(va[i].z);
    }
    mvalue.x *= inv_length;
    mvalue.y *= inv_length;
    mvalue.z *= inv_length;
    return mvalue;
  }
  else {
    rtErr("Data type " + std::string(typeid(T).name()) + " is not a suitable 3-tuple.",
          "trivariateMean");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> double3 trivariateMean(const std::vector<T> &va) {
  return trivariateMean(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double3 trivariateMean(const Hybrid<T> &va) {
  return trivariateMean(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double4 quadrivariateMean(const T* va, const size_t length) {
  if (isHpcVectorType<T>() && getHpcVectorTypeSize<T>() == 4) {
    const double inv_length = 1.0 / static_cast<double>(length);
    double4 mvalue = { 0.0, 0.0, 0.0, 0.0 };
    for (size_t i = 0; i < length; i++) {
      mvalue.x += static_cast<double>(va[i].x);
      mvalue.y += static_cast<double>(va[i].y);
      mvalue.z += static_cast<double>(va[i].z);
      mvalue.w += static_cast<double>(va[i].w);
    }
    mvalue.x *= inv_length;
    mvalue.y *= inv_length;
    mvalue.z *= inv_length;
    mvalue.w *= inv_length;
    return mvalue;
  }
  else {
    rtErr("Data type " + std::string(typeid(T).name()) + " is not a suitable 4-tuple.",
          "quadrivariateMean");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> double4 quadrivariateMean(const std::vector<T> &va) {
  return quadrivariateMean(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double4 quadrivariateMean(const Hybrid<T> &va) {
  return quadrivariateMean(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> T maxValue(const T* va, const size_t length) {
  if (isScalarType<T>() == false) {
    rtErr("Data type " + std::string(typeid(T).name()) + " is not suitable for computing a "
          "maximum value.", "maxValue");
  }
  if (length == 0) {
    rtErr("There is no data to scan for a maximum value.", "maxValue");
  }
  T max_value = va[0];
  for (size_t i = 1; i < length; i++) {
    max_value = std::max(max_value, va[i]);
  }
  return max_value;
}

//-------------------------------------------------------------------------------------------------
template <typename T> T maxValue(const std::vector<T> &va) {
  return maxValue(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> T maxValue(const Hybrid<T> &va) {
  return maxValue(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> T minValue(const T* va, const size_t length) {
  if (isScalarType<T>() == false) {
    rtErr("Data type " + std::string(typeid(T).name()) + " is not suitable for computing a "
          "minimum value.", "minValue");
  }
  if (length == 0) {
    rtErr("There is no data to scan for a minimum value.", "minValue");
  }
  T min_value = va[0];
  for (size_t i = 1; i < length; i++) {
    min_value = std::min(min_value, va[i]);
  }
  return min_value;
}

//-------------------------------------------------------------------------------------------------
template <typename T> T minValue(const std::vector<T> &va) {
  return minValue(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> T minValue(const Hybrid<T> &va) {
  return minValue(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> T maxAbsValue(const T* va, const size_t length) {
  if (isScalarType<T>() == false) {
    rtErr("Data type " + std::string(typeid(T).name()) + " is not suitable for computing a "
          "maximum absolute value.", "maxAbsValue");
  }
  if (length == 0) {
    rtErr("There is no data to scan for a minimum value.", "minValue");
  }
  T maxabs_value = (va[0] < static_cast<T>(0)) ? -va[0] : va[0];
  for (size_t i = 1; i < length; i++) {
    maxabs_value = std::max(maxabs_value, (va[i] < static_cast<T>(0)) ? -va[i] : va[i]);
  }
  return maxabs_value;
}

//-------------------------------------------------------------------------------------------------
template <typename T> T maxAbsValue(const std::vector<T> &va) {
  return maxAbsValue(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> T maxAbsValue(const Hybrid<T> &va) {
  return maxAbsValue(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double pearson(const T* va, const T* vb, const size_t length) {
  if (length < 2) {
    rtErr("There is insufficient data to correlate.", "minValue");
  }
  double sum_a  = 0.0;
  double sum_b  = 0.0;
  double sum_ab = 0.0;
  double sum_aa = 0.0;
  double sum_bb = 0.0;
  for (size_t i = 0; i < length; i++) {
    const double value_a = va[i];
    const double value_b = vb[i];
    sum_a += value_a;
    sum_b += value_b;
    sum_ab += value_a * value_b;
    sum_aa += value_a * value_a;
    sum_bb += value_b * value_b;
  }
  const double dlength = static_cast<double>(length);
  const double sq = (dlength*sum_aa - sum_a*sum_a) * (dlength*sum_bb - sum_b*sum_b);
  if (sq < constants::tiny) {
    return 0.0;
  }
  return (dlength*sum_ab - sum_a*sum_b) / sqrt(sq);
}

//-------------------------------------------------------------------------------------------------
template <typename T> double pearson(const std::vector<T> &va, const std::vector<T> &vb) {
  vectorComparisonCheck(va, vb, "pearson");
  return pearson(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double pearson(const Hybrid<T> &va, const Hybrid<T> &vb) {
  vectorComparisonCheck(va, vb, "pearson");
  return pearson(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> void addScalarToVector(T* va, const size_t length, const T inc) {
  for (size_t i = 0; i < length; i++) {
    va[i] += inc;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void addScalarToVector(std::vector<T> *va, const T inc) {
  addScalarToVector(va->data(), va->size(), inc);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void addScalarToVector(Hybrid<T> *va, const T inc) {
  addScalarToVector(va->data(), va->size(), inc);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void elementwiseMultiply(T* va, const size_t length, const T factor) {
  for (size_t i = 0; i < length; i++) {
    va[i] *= factor;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void elementwiseMultiply(std::vector<T> *va, const T factor) {
  elementwiseMultiply(va->data(), va->size(), factor);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void elementwiseMultiply(Hybrid<T> *va, const T factor) {
  elementwiseMultiply(va->data(), va->size(), factor);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void elementwiseDivide(T* va, const size_t length, const T factor) {
  for (size_t i = 0; i < length; i++) {
    va[i] /= factor;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void elementwiseDivide(std::vector<T> *va, const T factor) {
  elementwiseDivide(va->data(), va->size(), factor);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void elementwiseDivide(Hybrid<T> *va, const T factor) {
  elementwiseDivide(va->data(), va->size(), factor);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void crossProduct(const T* va, const T* vb, T* vc) {
  vc[0] = (va[1] * vb[2]) - (va[2] * vb[1]);
  vc[1] = (va[2] * vb[0]) - (va[0] * vb[2]);
  vc[2] = (va[0] * vb[1]) - (va[1] * vb[0]);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void crossProduct(const std::vector<T> &va, const std::vector<T> &vb,
                                        std::vector<T> *vc) {

  // Bounds check
  if (va.size() < 3 || vb.size() < 3) {
    rtErr("Vectors of size " + std::to_string(va.size()) + " and " + std::to_string(vb.size()) +
          "cannot produce a cross product.", "crossProduct");
  }
  else if (vc->size() < 3) {
    rtErr("A vector of size " + std::to_string(vc->size()) + " cannot hold a cross product.",
          "crossProduct");
  }

  // Set a pointer to the data of vector C and compute
  T* c_ptr = vc->data();
  c_ptr[0] = (va[1] * vb[2]) - (va[2] * vb[1]);
  c_ptr[1] = (va[2] * vb[0]) - (va[0] * vb[2]);
  c_ptr[2] = (va[0] * vb[1]) - (va[1] * vb[0]);
}

//-------------------------------------------------------------------------------------------------
template <typename T> T crossProduct(const T va, const T vb) {
  T vc;
  vc.x = (va.y * vb.z) - (va.z * vb.y);
  vc.y = (va.z * vb.x) - (va.x * vb.z);
  vc.z = (va.x * vb.y) - (va.y * vb.x);
  return vc;
}

//-------------------------------------------------------------------------------------------------
template <typename T> T hypotenuse(const T a, const T b) {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const T absa = std::abs(a);
  const T absb = std::abs(b);
  if (ct == double_type_index) {
    if (absa > absb) {
      return absa * sqrt(1.0 + ((absb / absa) * (absb / absa)));
    }
    else {
      return absb * sqrt(1.0 + ((absa / absb) * (absa / absb)));
    }
  }
  else {
    if (absa > absb) {
      return absa * sqrtf(1.0f + ((absb / absa) * (absb / absa)));
    }
    else {
      return absb * sqrtf(1.0f + ((absa / absb) * (absa / absb)));
    }
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
template <typename T> double magnitude(const T* va, const size_t length) {
  double result = 0.0;
  for (size_t i = 0; i < length; i++) {
    result += va[i] * va[i];
  }
  return sqrt(result);
}

//-------------------------------------------------------------------------------------------------
template <typename T> double magnitude(const std::vector<T> &va) {
  return magnitude(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double magnitude(const Hybrid<T> &va) {
  return magnitude(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double normalize(T* va, const size_t length) {
  const double magvec = magnitude(va, length);
  const double invmag = 1.0 / magvec;
  for (size_t i = 0; i < length; i++) {
    va[i] *= invmag;
  }
  return magvec;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double normalize(std::vector<T> *va) {
  return normalize(va->data(), va->size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double normalize(Hybrid<T> *va) {
  return normalize(va->data(), va->size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double dot(const T* va, const T* vb, const size_t length) {
  double result = 0.0;
  for (size_t i = 0; i < length; i++) {
    result += va[i] * vb[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double dot(const std::vector<T> &va, const std::vector<T> &vb) {
  vectorComparisonCheck(va, vb, "dot");
  return dot(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double dot(const Hybrid<T> &va, const Hybrid<T> &vb) {
  vectorComparisonCheck(va, vb, "dot");
  return dot(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double angleBetweenVectors(const T* va, const T* vb, const size_t length) {

  // Do this in slightly more optimized fashion
  double dotp_acc = 0.0;
  double maga_acc = 0.0;
  double magb_acc = 0.0;
  for (size_t i = 0; i < length; i++) {
    const double vai = va[i];
    const double vbi = vb[i];
    dotp_acc += vai * vbi;
    maga_acc += vai * vai;
    magb_acc += vbi * vbi;
  }
  const double mag2_ab = maga_acc * magb_acc;
  if (fabs(mag2_ab) < constants::tiny) {
    rtErr("One or both vectors are of close to zero length.", "angleBetweenVectors");
  }
  const double acos_arg = dotp_acc / (sqrt(maga_acc) * sqrt(magb_acc));
  if (acos_arg >= 1.0) {
    return 0.0;
  }
  else if (acos_arg <= -1.0) {
    return symbols::pi;
  }
  else {
    return acos(acos_arg);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> double angleBetweenVectors(const std::vector<T> &va,
                                                 const std::vector<T> &vb) {
  if (va.size() != vb.size()) {
    rtErr("Vectors of differing sizes, " + std::to_string(va.size()) + " and " +
          std::to_string(vb.size()) + ", cannot produce an angle value.", "angleBetweenVectors");
  }
  return angleBetweenVectors(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double angleBetweenVectors(const Hybrid<T> &va, const Hybrid<T> &vb) {
  if (va.size() != vb.size()) {
    rtErr("Hybrid objects " + std::string(va.getLabel().name) + " and " +
          std::string(vb.getLabel().name) + " have differing sizes, " + std::to_string(va.size()) +
          " and " + std::to_string(vb.size()) + ", and therefore cannot produce an angle value.",
          "angleBetweenVectors");
  }
  return angleBetweenVectors(va.data(), vb.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> void project(const T* va, const T* vb, T* vc, const size_t length) {
  const double mag_vb = magnitude(vb, length);
  const double dp_val = dot(va, vb, length) / (mag_vb * mag_vb);
  for (size_t i = 0; i < length; i++) {
    vc[i] = vb[i] * dp_val;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void project(const std::vector<T> &va, const std::vector<T> &vb,
                                   std::vector<T> *vc) {
  const size_t va_len = va.size();
  if (vb.size() != va_len || vc->size() != va_len) {
    rtErr("Vectors of length " + std::to_string(va_len) + " (vA), " + std::to_string(vb.size()) +
          " (vB), and " + std::to_string(vc->size()) + " (vC) are invalid for computing the "
          "projection of vA onto vB and storing it as vC.", "project");
  }
  project(va.data(), vb.data(), vc->data(), va_len);
}
  
//-------------------------------------------------------------------------------------------------
template <typename T> void project(const Hybrid<T> &va, const Hybrid<T> &vb, Hybrid<T> *vc) {
  const size_t va_len = va.size();
  if (vb.size() != va_len || vc->size() != va_len) {
    rtErr("Vectors of length " + std::to_string(va_len) + " (vA), " + std::to_string(vb.size()) +
          " (vB), and " + std::to_string(vc->size()) + " (vC) are invalid for computing the "
          "projection of vA (" + std::string(va.getLabel().name) + ") onto vB (" +
          std::string(vb.getLabel().name) + ") and storing it as vC (" +
          std::string(vc->getLabel().name) + ").", "project");
  }
  project(va.data(), vb.data(), vc->data(), va_len);
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> project(const std::vector<T> &va, const std::vector<T> &vb) {
  vectorComparisonCheck(va, vb, "project");
  const size_t va_len = va.size();
  const double mag_vb = magnitude(vb.data(), va_len);
  const double dp_val = dot(va.data(), vb.data(), va_len) / (mag_vb * mag_vb);
  std::vector<T> result(va_len);
  for (size_t i = 0; i < va_len; i++) {
    result[i] = vb[i] * dp_val;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> void perpendicularComponent(const T* va, const T* vb, T* result,
                                                  const size_t length) {
  project(va, vb, result, length);
  for (size_t i = 0; i < length; i++) {
    result[i] = va[i] - result[i];
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void perpendicularComponent(const std::vector<T> &va,
                                                  const std::vector<T> &vb,
                                                  std::vector<T> *result) {
  vectorComparisonCheck(va, vb, "perpendicularComponent");
  vectorComparisonCheck(va, result, "perpendicularComponent");
  perpendicularComponent(va.data(), vb.data(), result->data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> perpendicularComponent(const std::vector<T> &va,
                                                            const std::vector<T> &vb) {
  vectorComparisonCheck(va, vb, "perpendicularComponent");
  std::vector<double> result = project(va, vb);
  const size_t va_dim = va.size();
  for (size_t i = 0; i < va_dim; i++) {
    result[i] = va[i] - result[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double pointPlaneDistance(const T* va, const T* vb, const T* pt_pos) {
  double unit_normal[3], pt_plane_displacement[3];
  crossProduct(va, vb, unit_normal);
  normalize(unit_normal, 3);
  project(unit_normal, pt_pos, pt_plane_displacement, 3);
  return magnitude(pt_plane_displacement, 3);
}

//-------------------------------------------------------------------------------------------------
template <typename T> double pointPlaneDistance(const std::vector<T> &va, const std::vector<T> &vb,
                                                const std::vector<T> &pt_pos) {
  if (va.size() != 3LLU || vb.size() != 3LLU || pt_pos.size() != 3LLU) {
    rtErr("The distance from a point to a plane is computed in three dimensions.  Dimensions of "
          "vectors provided = [ " + std::to_string(va.size()) + ", " + std::to_string(vb.size()) +
          ", " + std::to_string(pt_pos.size()) + " ].", "pointPlaneDistance");
  }
  return pointPlaneDistance(va.data(), vb.data(), pt_pos.data());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double pointPlaneDistance(const Hybrid<T> &va, const Hybrid<T> &vb,
                                                const Hybrid<T> &pt_pos) {
  if (va.size() != 3LLU || vb.size() != 3LLU || pt_pos.size() != 3LLU) {
    rtErr("The distance from a point to a plane is computed in three dimensions.  Dimensions of "
          "Hybrid objects " + std::string(va.getLabel().name) + " and " +
          std::string(vb.getLabel().name) + " are [ " + std::to_string(va.size()) + ", " +
          std::to_string(vb.size()) + ", " + std::to_string(pt_pos.size()) + " ].",
          "pointPlaneDistance");
  }
  return pointPlaneDistance(va.data(), vb.data(), pt_pos.data());
}

//-------------------------------------------------------------------------------------------------
template <typename T> void reportBinLimitError(const std::string &desc, const T value,
                                               const ExceptionResponse policy) {
  switch (policy) {
  case ExceptionResponse::DIE:
    if (isSignedIntegralScalarType<T>()) {
      rtErr("A value of " + std::to_string(static_cast<llint>(value)) + 
            " is off the " + desc + " end of the range.", "findBin");
    }
    else if (isUnsignedIntegralScalarType<T>()) {
      rtErr("A value of " + std::to_string(static_cast<ullint>(value)) + 
            " is off the " + desc + " end of the range.", "findBin");
    }
    else if (isFloatingPointScalarType<T>()) {
      rtErr("A value of " +
            realToString(static_cast<double>(value), 11, 4, NumberFormat::SCIENTIFIC) +
            " is off the " + desc + " end of the range.", "findBin");
    }
    break;
  case ExceptionResponse::WARN:
    if (isSignedIntegralScalarType<T>()) {
      rtWarn("A value of " + std::to_string(static_cast<llint>(value)) + 
             " is off the " + desc + " end of the range.", "findBin");
    }
    else if (isUnsignedIntegralScalarType<T>()) {
      rtWarn("A value of " + std::to_string(static_cast<ullint>(value)) + 
             " is off the " + desc + " end of the range.", "findBin");
    }
    else if (isFloatingPointScalarType<T>()) {
      rtWarn("A value of " +
             realToString(static_cast<double>(value), 11, 4, NumberFormat::SCIENTIFIC) +
             " is off the " + desc + " end of the range.", "findBin");
    }
    break;
  case ExceptionResponse::SILENT:
    break;
  }
}
  
//-------------------------------------------------------------------------------------------------
template <typename T> int findBin(const T* limits, const T value, const int length,
                                  const ExceptionResponse policy) {
  if (length == 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A zero-length set of bin limits was supplied.", "findBin");
    case ExceptionResponse::WARN:
      rtWarn("A zero-length set of bin limits was supplied.", "findBin");
      return -1;
    case ExceptionResponse::SILENT:
      return -1;
    }
  }
  if (value < limits[0]) {
    reportBinLimitError("left", value, policy);
    return -1;
  }
  if (value >= limits[length]) {
    reportBinLimitError("right", value, policy);
    return length + 1;
  }
  int lguess = 0;
  int hguess = length;
  while (lguess < hguess - 1LLU) {

    // Choose a residue intermediate between the lower and upper bounds
    const int mguess = lguess + ((hguess - lguess) / 2);
    if (value >= limits[mguess + 1LLU]) {
      lguess = mguess;
    }
    else if (value < limits[mguess]) {
      hguess = mguess;
    }
    else {
      return mguess;
    }
  }
  return lguess;
}

//-------------------------------------------------------------------------------------------------
template <typename T> int findBin(const std::vector<T> &limits, const T value,
                                  const ExceptionResponse policy) {
  return findBin(limits.data(), value, limits.size() - 1LLU, policy);
}

//-------------------------------------------------------------------------------------------------
template <typename T> int findBin(const Hybrid<T> &limits, const T value,
                                  const ExceptionResponse policy) {
  return findBin(limits.data(), value, limits.size() - 1LLU, policy);
}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t locateValue(const T* vdata, const T value, const size_t length,
                                         const DataOrder format) {
  if (length == 0LLU) {
    rtErr("Unable to search an array of zero length.", "locateValue");
  }
  switch (format) {
  case DataOrder::ASCENDING:
    if (value >= vdata[0] && value <= vdata[length - 1]) {
      size_t min_pos = 0LLU;
      size_t max_pos = length;
      while (max_pos - min_pos > 1) {
        size_t mid_pos = (max_pos + min_pos) / 2;
        if (vdata[mid_pos] > value) {
          max_pos = mid_pos;
        }
        else {
          min_pos = mid_pos;
        }
      }
      return (vdata[min_pos] == value) ? min_pos : length;
    }
    break;
  case DataOrder::DESCENDING:
    if (value <= vdata[0] && value >= vdata[length - 1]) {
      size_t min_pos = 0LLU;
      size_t max_pos = length;
      while (max_pos - min_pos > 1) {
        size_t mid_pos = (max_pos + min_pos) / 2;
        if (vdata[mid_pos] < value) {
          max_pos = mid_pos;
        }
        else {
          min_pos = mid_pos;
        }
      }
      return (vdata[min_pos] == value) ? min_pos : length;
    }
    break;
  case DataOrder::NONE:
    for (size_t i = 0; i < length; i++) {
      if (vdata[i] == value) {
        return i;
      }
    }
    break;
  }
  return length;
}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t locateValue(const std::vector<T> &vdata, const T value,
                                         const DataOrder format) {
  return locateValue(vdata.data(), value, vdata.size(), format);
}
  
//-------------------------------------------------------------------------------------------------
template <typename T> size_t locateValue(const Hybrid<T> &vdata, const T value,
                                         const DataOrder format) {
  return locateValue(vdata.data(), value, vdata.size(), format);
}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t locateValue(const Approx &value, const T* vdata, const size_t length,
                                         const DataOrder format) {
  if (length == 0LLU) {
    rtErr("Unable to search an array of zero length.", "locateValue");
  }
  switch (format) {
  case DataOrder::ASCENDING:
    if (value >= vdata[0] && value <= vdata[length - 1]) {
      size_t min_pos = 0LLU;
      size_t max_pos = length;
      while (max_pos - min_pos > 1) {
        size_t mid_pos = (max_pos + min_pos) / 2;
        if (vdata[mid_pos] > value) {
          max_pos = mid_pos;
        }
        else {
          min_pos = mid_pos;
        }
      }
      return (vdata[min_pos] == value) ? min_pos : length;
    }
    break;
  case DataOrder::DESCENDING:
    if (value <= vdata[0] && value >= vdata[length - 1]) {
      size_t min_pos = 0LLU;
      size_t max_pos = length;
      while (max_pos - min_pos > 1) {
        size_t mid_pos = (max_pos + min_pos) / 2;
        if (vdata[mid_pos] < value) {
          max_pos = mid_pos;
        }
        else {
          min_pos = mid_pos;
        }
      }
      return (vdata[min_pos] == value) ? min_pos : length;
    }
    break;
  case DataOrder::NONE:
    for (size_t i = 0; i < length; i++) {
      if (vdata[i] == value) {
        return i;
      }
    }
    break;
  }
  return length;
}
  
//-------------------------------------------------------------------------------------------------
template <typename T> size_t locateValue(const Approx &value, const std::vector<T> &vdata,
                                         const DataOrder format) {
  return locateValue(value, vdata.data(), vdata.size(), format);
}
  
//-------------------------------------------------------------------------------------------------
template <typename T> size_t locateValue(const Approx &value, const Hybrid<T> &vdata,
                                         const DataOrder format) {
  return locateValue(value, vdata.data(), vdata.size(), format);
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> reduceUniqueValues(const std::vector<T> &va) {
  std::vector<T> result(va);
  reduceUniqueValues(&result);
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> void reduceUniqueValues(std::vector<T> *va) {

  // Return immediately if there is no data
  if (va->size() == 0) {
    return;
  }

  // If the type is integral, test whether the data might be well clustered
  T* va_ptr = va->data();
  const size_t nval = va->size();
  if (isSignedIntegralScalarType<T>() || isUnsignedIntegralScalarType<T>()) {
    T va_min = va_ptr[0];
    T va_max = va_ptr[0];
    for (size_t i = 1LLU; i < nval; i++) {
      va_min = std::min(va_min, va_ptr[i]);
      va_max = std::max(va_max, va_ptr[i]);
    }
    const size_t max_unique = va_max - va_min + 1LLU;
    if (max_unique < nval) {
      std::vector<bool> presence(max_unique, false);
      for (size_t i = 0LLU; i < nval; i++) {
        presence[va_ptr[i] - va_min] = true;
      }
      size_t j = 0LLU;
      T ti = 0;
      for (size_t i = 0LLU; i < max_unique; i++) {
        if (presence[i]) {
          va_ptr[j] = va_min + ti;
          j++;
        }
        ti++;
      }
      va->resize(j);
      va->shrink_to_fit();      
      return;
    }
  }

  // Take the standard sorting approach
  std::sort(va->begin(), va->end(), [](T a, T b) { return a < b; });
  size_t nunique = 0LLU;
  T last_unique;
  for (size_t i = 0LLU; i < nval; i++) {
    if (i == 0LLU || va_ptr[i] != last_unique) {
      last_unique = va_ptr[i];
      va_ptr[nunique] = last_unique;
      nunique++;
    }
  }
  va->resize(nunique);
  va->shrink_to_fit();
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<ValueWithCounter<T>>
findUnmatchedValues(const std::vector<T> &va, const std::vector<T> &vb,
                    const UniqueValueHandling check_repeats) {
  const size_t length_a = va.size();
  const size_t length_b = vb.size();
  std::vector<T> va_duplicate(va);
  std::vector<T> vb_duplicate(vb);
  std::sort(va_duplicate.begin(), va_duplicate.end(), [](T a, T b) { return a < b; });
  std::sort(vb_duplicate.begin(), vb_duplicate.end(), [](T a, T b) { return a < b; });
  std::vector<ValueWithCounter<T>> unique_a;
  std::vector<ValueWithCounter<T>> unique_b;
  for (size_t i = 0; i < length_a; i++) {
    if (i == 0 || va_duplicate[i - 1] != va_duplicate[i]) {
      unique_a.push_back({va_duplicate[i], 1});
    }
    else {
      unique_a[unique_a.size() - 1LLU].count += 1;
    }
  }
  for (size_t i = 0; i < length_b; i++) {
    if (i == 0 || vb_duplicate[i - 1] != vb_duplicate[i]) {
      unique_b.push_back({vb_duplicate[i], 1});
    }
    else {
      unique_b[unique_b.size() - 1LLU].count += 1;
    }
  }
  const size_t nua = unique_a.size();
  const size_t nub = unique_b.size();
  size_t min_j = 0;
  std::vector<ValueWithCounter<T>> result;
  const bool counts_matter = (check_repeats == UniqueValueHandling::CONFIRM_ALL_COPIES);
  for (size_t i = 0; i < nua; i++) {
    const T aval = unique_a[i].value;
    const int acount = unique_a[i].count;
    bool found = false;
    for (size_t j = min_j; j < nub; j++) {
      if (unique_b[j].value == aval) {
        min_j += (j == min_j);
        if (counts_matter && unique_b[j].count != acount) {
          result.push_back({aval, acount - unique_b[j].count});
        }
        found = true;
      }
    }
    if (! found) {
      result.push_back({aval, acount});
    }
  }
  min_j = 0;
  for (size_t i = 0; i < nub; i++) {
    const T bval = unique_b[i].value;
    const int bcount = unique_b[i].count;
    bool found = false;
    for (size_t j = min_j; j < nua; j++) {
      if (unique_a[j].value == bval) {
        min_j += (j == min_j);
        found = true;
      }
    }
    if (! found) {
      result.push_back({bval, -bcount});
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<ValueWithCounter<T>>
findUnmatchedValues(const T* va, const T* vb, const size_t length_a, const size_t length_b,
                    const UniqueValueHandling check_repeats) {
  std::vector<T> tva(length_a);
  std::vector<T> tvb(length_b);
  return findUnmatchedValues(tva, tvb, check_repeats);
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<ValueWithCounter<T>>
findUnmatchedValues(const Hybrid<T> &va, const Hybrid<T> &vb,
                    const UniqueValueHandling check_repeats) {
  std::vector<T> tva = va.readHost();
  std::vector<T> tvb = vb.readHost();
  return findUnmatchedValues(tva, tvb, check_repeats);
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<bool> colorVectorMask(const std::vector<T> &mask) {
  const size_t nmask = mask.size();

  // Determine the size of the result based on the largest entry of the mask
  T nrslt = 0;
  for (size_t i = 0; i < mask[i]; i++) {
    if (mask[i] < 0) {
      rtErr("A negative index (" + std::to_string(mask[i]) + ") cannot be marked on a vector "
            "mask.", "colorVectorMask");
    }
    nrslt = std::max(nrslt, mask[i]);
  }
  std::vector<bool> result(nrslt, false);
  for (size_t i = 0; i < mask[i]; i++) {
    result[mask[i]] = true;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<bool> colorVectorMask(const std::vector<T> &mask, const size_t length) {

  // Check that the largest entry of the mask does not exceed the specified length when coloring
  // the result.
  std::vector<bool> result(length, false);
  for (size_t i = 0; i < mask[i]; i++) {
    if (mask[i] < 0) {
      rtErr("A negative index (" + std::to_string(mask[i]) + ") cannot be marked on a vector "
            "mask.", "colorVectorMask");
    }
    if (mask[i] >= length) {
      rtErr("An array of specified length " + std::to_string(length) + " cannot mark an index " +
            "of " + std::to_string(mask[i]) + ".", "colorVectorMask");
    }
    result[mask[i]] = true;
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
template <typename T> void colorVectorMask(std::vector<bool> *result, const std::vector<T> &mask) {
  const size_t nmask = mask.size();
  const size_t nrslt = result->size();
  for (size_t i = 0; i < nmask; i++) {
    if (mask[i] < 0 || mask[i] >= nrslt) {
      rtErr("An array of " + std::to_string(nrslt) + " cannot mark an index of " +
            std::to_string(mask[i]) + ".", "colorVectorMask");
    }
    result->at(mask[i]) = true;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> tileVector(const std::vector<T> &va, const size_t nrep) {
  std::vector<T> result;
  const size_t n_va = va.size();
  result.reserve(n_va * nrep);
  for (size_t pos = 0LLU; pos < nrep; pos++) {	
    for (size_t i = 0LLU; i < n_va; i++) {
      result.emplace_back(va[i]);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> tileVector(const std::vector<T> &va, const std::vector<int> &tidx,
                          const size_t nrep) {
  std::vector<T> result;
  const size_t n_tidx = tidx.size();
  result.reserve(n_tidx * nrep);
  for (size_t pos = 0LLU; pos < nrep; pos++) {
    for (size_t i = 0LLU; i < n_tidx; i++) {
      result.emplace_back(va[tidx[i]]);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> tileVector(const std::vector<T> &va, const std::vector<size_t> &tidx,
                          const size_t nrep) {
  std::vector<T> result;
  const size_t n_tidx = tidx.size();
  result.reserve(n_tidx * nrep);
  for (size_t pos = 0LLU; pos < nrep; pos++) {
    for (size_t i = 0LLU; i < n_tidx; i++) {
      result.emplace_back(va[tidx[i]]);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tout>
std::vector<Tout> convertData(std::vector<int95_t> *input, const double input_scale,
                              const double output_scale) {
  const int95_t* inp_ptr = input->data();
  const size_t npt = input->size();
  std::vector<Tout> output(npt);
  if (isFloatingPointScalarType<Tout>()) {
    const double conv_fac = output_scale / input_scale;
    for (size_t i = 0; i < npt; i++) {
      output[i] = hostInt95ToDouble(inp_ptr[i]) * conv_fac;
    }
  }
  else if (isSignedIntegralScalarType<Tout>() || isUnsignedIntegralScalarType<Tout>()) {
    const size_t ct_out = std::type_index(typeid(Tout)).hash_code();
    if (ct_out == llint_type_index || ct_out == ullint_type_index) {
      const int input_scale_bits  = log2(input_scale);
      const int output_scale_bits = log2(output_scale);
      for (size_t i = 0; i < npt; i++) {
        const int95_t adj = hostChangeFPBits(inp_ptr[i], input_scale_bits, output_scale_bits);

        // Only the 64-bit part can be copied into the output.  Conversions between int95_t
        // vectors should be handled by running hostChangeFPBits() over the one vector.
        output[i] = adj.x;
      }
    }
    else {

      // Double-precision can handle other conversions without loss of information.
      const double conv_fac = output_scale / input_scale;
      for (size_t i = 0; i < npt; i++) {
        output[i] = llround(hostInt95ToDouble(inp_ptr[i]) * conv_fac);
      }
    }
  }
  input->resize(0);
  input->shrink_to_fit();
  return output;
}

//-------------------------------------------------------------------------------------------------
template <typename Tout>
std::vector<Tout> convertData(std::vector<llint> *primary, std::vector<int> *overflow,
                              const double input_scale, const double output_scale) {
  const llint* prim_ptr = primary->data();
  const int* ovrf_ptr = overflow->data();
  const size_t npt = primary->size();
  if (overflow->size() != npt) {
    rtErr("Input primary and overflow arrays must be of equal size (provided " +
          std::to_string(npt) + " and " + std::to_string(overflow->size()) + ").", "convertData");
  }
  std::vector<Tout> output(npt);
  if (isFloatingPointScalarType<Tout>()) {
    const double conv_fac = output_scale / input_scale;
    for (size_t i = 0; i < npt; i++) {
      output[i] = hostInt95ToDouble(prim_ptr[i], ovrf_ptr[i]) * conv_fac;
    }
  }
  else if (isSignedIntegralScalarType<Tout>() || isUnsignedIntegralScalarType<Tout>()) {
    const size_t ct_out = std::type_index(typeid(Tout)).hash_code();
    if (ct_out == llint_type_index || ct_out == ullint_type_index) {
      const int input_scale_bits  = log2(input_scale);
      const int output_scale_bits = log2(output_scale);
      for (size_t i = 0; i < npt; i++) {
        const int95_t orig = { prim_ptr[i], ovrf_ptr[i] };
        const int95_t adj = hostChangeFPBits(orig, input_scale_bits, output_scale_bits);

        // Only the 64-bit part can be copied into the output.  Conversions between int95_t
        // vectors should be handled by running hostChangeFPBits() over the one vector.
        output[i] = adj.x;
      }
    }
    else {

      // Double-precision can handle other conversions without loss of information.
      const double conv_fac = output_scale / input_scale;
      for (size_t i = 0; i < npt; i++) {
        output[i] = llround(hostInt95ToDouble(prim_ptr[i], ovrf_ptr[i]) * conv_fac);
      }
    }
  }
  primary->resize(0);
  primary->shrink_to_fit();
  overflow->resize(0);
  overflow->shrink_to_fit();
  return output;
}

//-------------------------------------------------------------------------------------------------
template <typename Tin>
std::vector<int95_t> convertData(std::vector<Tin> *input, const double input_scale,
                                 const double output_scale) {
  const Tin* inp_ptr = input->data();
  const size_t npt = input->size();
  std::vector<int95_t> output(npt);
  if (isFloatingPointScalarType<Tin>()) {
    const double conv_fac = output_scale / input_scale;
    for (size_t i = 0; i < npt; i++) {
      output[i] = hostDoubleToInt95(input[i] * conv_fac);
    }
  }
  else if (isSignedIntegralScalarType<Tin>() || isUnsignedIntegralScalarType<Tin>()) {
    const size_t ct_in = std::type_index(typeid(Tin)).hash_code();
    if (ct_in == llint_type_index || ct_in == ullint_type_index) {
      const int input_scale_bits  = log2(input_scale);
      const int output_scale_bits = log2(output_scale);
      for (size_t i = 0; i < npt; i++) {
        const int95_t orig = { inp_ptr[i], 0 };
        const int95_t adj = hostChangeFPBits(orig, input_scale_bits, output_scale_bits);
        output[i] = adj;
      }
    }
    else {
      const double conv_fac = output_scale / input_scale;
      for (size_t i = 0; i < npt; i++) {
        output[i] = hostDoubleToInt95(static_cast<double>(inp_ptr[i]) * conv_fac);
      }
    }
  }
  input->resize(0);
  input->shrink_to_fit();
  return output;
}

//-------------------------------------------------------------------------------------------------
template <typename Tin, typename Tout>
std::vector<Tout> convertData(std::vector<Tin> *input, const double input_scale,
                              const double output_scale) {
  const Tin* inp_ptr = input->data();
  const size_t npt = input->size();
  std::vector<Tout> output(npt);
  const size_t ct_in  = std::type_index(typeid(Tin)).hash_code();
  const size_t ct_out = std::type_index(typeid(Tout)).hash_code();
  if ((ct_in  == llint_type_index || ct_in  == ullint_type_index) &&
      (ct_out == llint_type_index || ct_out == ullint_type_index)) {
    if (input_scale > output_scale) {
      const int bit_delta = round(log2(input_scale / output_scale));
      for (size_t i = 0; i < npt; i++) {
        output[i] = (inp_ptr[i] >> bit_delta);
      }
    }
    else if (output_scale > input_scale) {
      const int bit_delta = round(log2(output_scale / input_scale));
      for (size_t i = 0; i < npt; i++) {
        output[i] = (inp_ptr[i] << bit_delta);
      }
    }
    else {
      for (size_t i = 0; i < npt; i++) {
        output[i] = inp_ptr[i];
      }
    }
  }
  else {
    const double conv_fac = output_scale / input_scale;
    if (isFloatingPointScalarType<Tout>()) {
      for (size_t i = 0; i < npt; i++) {
        output[i] = static_cast<double>(inp_ptr[i]) * conv_fac;
      }
    }
    else {
      for (size_t i = 0; i < npt; i++) {
        output[i] = llround(static_cast<double>(inp_ptr[i]) * conv_fac);
      }
    }
  }
  input->resize(0);
  input->shrink_to_fit();
  return output;
}

//-------------------------------------------------------------------------------------------------
template <typename Tout>
Hybrid<Tout> convertData(Hybrid<int95_t> *input, const double input_scale,
                         const double output_scale) {
  const int95_t* inp_ptr = input->data();
  const size_t npt = input->size();
  Hybrid<Tout> output(npt, input->getLabel().name, input->getFormat());
  Tout* out_ptr = output.data();
  if (isFloatingPointScalarType<Tout>()) {
    const double conv_fac = output_scale / input_scale;
    for (size_t i = 0; i < npt; i++) {
      out_ptr[i] = hostInt95ToDouble(inp_ptr[i]) * conv_fac;
    }
  }
  else if (isSignedIntegralScalarType<Tout>() || isUnsignedIntegralScalarType<Tout>()) {
    const size_t ct_out = std::type_index(typeid(Tout)).hash_code();
    if (ct_out == llint_type_index || ct_out == ullint_type_index) {
      const int input_scale_bits  = log2(input_scale);
      const int output_scale_bits = log2(output_scale);
      for (size_t i = 0; i < npt; i++) {
        const int95_t adj = hostChangeFPBits(inp_ptr[i], input_scale_bits, output_scale_bits);

        // Only the 64-bit part can be copied into the output.  Conversions between int95_t
        // vectors should be handled by running hostChangeFPBits() over the one vector.
        out_ptr[i] = adj.x;
      }
    }
    else {

      // Double-precision can handle other conversions without loss of information.
      const double conv_fac = output_scale / input_scale;
      for (size_t i = 0; i < npt; i++) {
        out_ptr[i] = llround(hostInt95ToDouble(inp_ptr[i]) * conv_fac);
      }
    }
  }
  input->resize(0);
  input->shrinkToFit();
  return output;
}

//-------------------------------------------------------------------------------------------------
template <typename Tout>
Hybrid<Tout> convertData(Hybrid<llint> *primary, Hybrid<int> *overflow,
                         const double input_scale, const double output_scale) {
  const llint* prim_ptr = primary->data();
  const int* ovrf_ptr = overflow->data();
  const size_t npt = primary->size();
  if (overflow->size() != npt) {
    rtErr("Input primary and overflow arrays must be of equal size (provided " +
          std::to_string(npt) + " and " + std::to_string(overflow->size()) + ").", "convertData");
  }
  Hybrid<Tout> output(npt, primary->getLabel().name, primary->getFormat());
  Tout* out_ptr = output.data();
  if (isFloatingPointScalarType<Tout>()) {
    const double conv_fac = output_scale / input_scale;
    for (size_t i = 0; i < npt; i++) {
      out_ptr[i] = hostInt95ToDouble(prim_ptr[i], ovrf_ptr[i]) * conv_fac;
    }
  }
  else if (isSignedIntegralScalarType<Tout>() || isUnsignedIntegralScalarType<Tout>()) {
    const size_t ct_out = std::type_index(typeid(Tout)).hash_code();
    if (ct_out == llint_type_index || ct_out == ullint_type_index) {
      const int input_scale_bits  = log2(input_scale);
      const int output_scale_bits = log2(output_scale);
      for (size_t i = 0; i < npt; i++) {
        const int95_t orig = { prim_ptr[i], ovrf_ptr[i] };
        const int95_t adj = hostChangeFPBits(orig, input_scale_bits, output_scale_bits);

        // Only the 64-bit part can be copied into the output.  Conversions between int95_t
        // vectors should be handled by running hostChangeFPBits() over the one vector.
        out_ptr[i] = adj.x;
      }
    }
    else {

      // Double-precision can handle other conversions without loss of information.
      const double conv_fac = output_scale / input_scale;
      for (size_t i = 0; i < npt; i++) {
        out_ptr[i] = llround(hostInt95ToDouble(prim_ptr[i], ovrf_ptr[i]) * conv_fac);
      }
    }
  }
  primary->resize(0);
  primary->shrinkToFit();
  overflow->resize(0);
  overflow->shrinkToFit();
  return output;
}

//-------------------------------------------------------------------------------------------------
template <typename Tin>
Hybrid<int95_t> convertData(Hybrid<Tin> *input, const double input_scale,
                            const double output_scale) {
  const Tin* inp_ptr = input->data();
  const size_t npt = input->size();
  Hybrid<int95_t> output(npt, input->getLabel().name, input->getFormat());
  int95_t* out_ptr = output.data();
  if (isFloatingPointScalarType<Tin>()) {
    const double conv_fac = output_scale / input_scale;
    for (size_t i = 0; i < npt; i++) {
      out_ptr[i] = hostDoubleToInt95(input[i] * conv_fac);
    }
  }
  else if (isSignedIntegralScalarType<Tin>() || isUnsignedIntegralScalarType<Tin>()) {
    const size_t ct_in = std::type_index(typeid(Tin)).hash_code();
    if (ct_in == llint_type_index || ct_in == ullint_type_index) {
      const int input_scale_bits  = log2(input_scale);
      const int output_scale_bits = log2(output_scale);
      for (size_t i = 0; i < npt; i++) {
        const int95_t orig = { inp_ptr[i], 0 };
        const int95_t adj = hostChangeFPBits(orig, input_scale_bits, output_scale_bits);
        out_ptr[i] = adj;
      }
    }
    else {
      const double conv_fac = output_scale / input_scale;
      for (size_t i = 0; i < npt; i++) {
        out_ptr[i] = hostDoubleToInt95(static_cast<double>(inp_ptr[i]) * conv_fac);
      }
    }
  }
  input->resize(0);
  input->shrinkToFit();
  return output;
}

//-------------------------------------------------------------------------------------------------
template <typename Tin, typename Tout>
Hybrid<Tout> convertData(Hybrid<Tin> *input, const double input_scale, const double output_scale) {
  const Tin* inp_ptr = input->data();
  const size_t npt = input->size();
  Hybrid<Tout> output(npt, input->getLabel().name, input->getFormat());
  Tout* out_ptr = output.data();
  const size_t ct_in  = std::type_index(typeid(Tin)).hash_code();
  const size_t ct_out = std::type_index(typeid(Tout)).hash_code();
  if ((ct_in  == llint_type_index || ct_in  == ullint_type_index) &&
      (ct_out == llint_type_index || ct_out == ullint_type_index)) {
    if (input_scale > output_scale) {
      const int bit_delta = round(log2(input_scale / output_scale));
      for (size_t i = 0; i < npt; i++) {
        out_ptr[i] = (inp_ptr[i] >> bit_delta);
      }
    }
    else if (output_scale > input_scale) {
      const int bit_delta = round(log2(output_scale / input_scale));
      for (size_t i = 0; i < npt; i++) {
        out_ptr[i] = (inp_ptr[i] << bit_delta);
      }
    }
    else {
      for (size_t i = 0; i < npt; i++) {
        out_ptr[i] = inp_ptr[i];
      }
    }
  }
  else {
    const double conv_fac = output_scale / input_scale;
    if (isFloatingPointScalarType<Tout>()) {
      for (size_t i = 0; i < npt; i++) {
        out_ptr[i] = static_cast<double>(inp_ptr[i]) * conv_fac;
      }
    }
    else {
      for (size_t i = 0; i < npt; i++) {
        out_ptr[i] = llround(static_cast<double>(inp_ptr[i]) * conv_fac);
      }
    }
  }
  input->resize(0);
  input->shrinkToFit();
  return output;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> applyAssociatedSort(const std::vector<T> &unsorted_data,
                                   const std::vector<int> &sorting_pattern) {
  return tileVector<T>(unsorted_data, sorting_pattern, 1);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> applyAssociatedSort(const std::vector<T> &unsorted_data,
                                   const std::vector<size_t> &sorting_pattern) {
  return tileVector<T>(unsorted_data, sorting_pattern, 1);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2>
std::vector<T> applyAssociatedSort(const std::vector<T> &unsorted_data,
                                   const std::vector<T2> &sorting_pattern) {
  const size_t plen = sorting_pattern.size();
  std::vector<T> result(plen);
  for (size_t i = 0; i < plen; i++) {
    result[i] = unsorted_data[sorting_pattern[i].y];
  }
  return result;
}

} // namespace stmath
} // namespace stormm
