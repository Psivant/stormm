// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace testing {

//-------------------------------------------------------------------------------------------------
template <typename T> Approx::Approx(const std::vector<T> &values_in, ComparisonType style_in,
                                     double tol_in) :
    values{},
    style{style_in},
    dtol{tol_in}
{
  // Check that the vector is of an acceptable type.  The vector will get converted to a
  // double-precision floating point representation, irrespective of its input type.  Very
  // large integers, as large as 53 bits, can be represented exactly, but the last eleven bits
  // are going to get approximated away.  A 64-bit signed or unsigned integer's low 11 or 12
  // bits will be lost if the number is exceptionally large.
  if (isScalarType<T>() == false) {
    rtErr("Data of type " + std::string(typeid(T).name()) +
          " cannot be processed for approximate, real-valued comparisons.", "Approx");
  }
  values = std::vector<double>(values_in.begin(), values_in.end());

  // Check that the tolerance is valid
  if (dtol < 0.0) {
    rtErr("A negative tolerance (" + std::to_string(dtol) + ") is nonsensical.", "Approx");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool verifyVectorApproxCompatibility(const std::vector<T> &test_values,
                                                           const Approx &cr) {

  // Check that the vector is of an acceptable type
  if (isScalarType<T>() == false) {
    rtErr("Real-valued equality comparisons for a vector of type " +
          std::string(typeid(T).name()) + " are not permitted.");
  }

  // Check that vectors are populated and their sizes match
  const size_t nval = cr.size();
  const size_t ntval = test_values.size();
  return (nval == ntval && nval > 0);
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool Approx::test(const std::vector<T> &test_values) const {
  if (verifyVectorApproxCompatibility(test_values, *this) == false) {
    return false;
  }
  const size_t nval = values.size();

  // Promote the test values to double precision to match the reference data type
  const std::vector<double> converted_values(test_values.begin(), test_values.end());

  // Declare the "accumulated deviation" outside the switch to avoid the need for extra scoping
  double acc_dev = 0.0;
  switch (style) {
  case ComparisonType::ABSOLUTE:
    for (size_t i = 0; i < nval; i++) {
      acc_dev = std::max(fabs(converted_values[i] - values[i]), acc_dev);
    }
    break;
  case ComparisonType::MEAN_UNSIGNED_ERROR:
    for (size_t i = 0; i < nval; i++) {
      acc_dev += fabs(converted_values[i] - values[i]);
    }
    acc_dev /= static_cast<double>(nval);
    break;
  case ComparisonType::RELATIVE:
    for (size_t i = 0; i < nval; i++) {
      if (std::abs(values[i]) > constants::tiny) {
        acc_dev += std::max(fabs((converted_values[i] - values[i]) / values[i]), acc_dev);
      }
      else {
        acc_dev += std::max(fabs((converted_values[i] - values[i]) / constants::tiny), acc_dev);
      }
    }
    break;
  case ComparisonType::RELATIVE_RMS_ERROR:
    {
      double abs_mean = 0.0;
      for (size_t i = 0; i < nval; i++) {
        abs_mean += fabs(values[i]);
        acc_dev += (converted_values[i] - values[i]) * (converted_values[i] - values[i]);
      }
      abs_mean /= static_cast<double>(nval);
      if (std::abs(abs_mean) > constants::tiny) {
        acc_dev = sqrt(acc_dev / static_cast<double>(nval)) / abs_mean;
      }
      else {
        acc_dev = sqrt(acc_dev / static_cast<double>(nval)) / constants::tiny;
      }
    }
    break;
  }
  return (acc_dev < dtol);
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool operator==(const std::vector<T> &tvec, const Approx &cr) {
  return cr.test(tvec);
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool operator==(const Approx &cr, const std::vector<T> &tvec) {
  return cr.test(tvec);
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool operator!=(const std::vector<T> &tvec, const Approx &cr) {
  return (cr.test(tvec) == false);
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool operator!=(const Approx &cr, const std::vector<T> &tvec) {
  return (cr.test(tvec) == false);
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool operator>(const std::vector<T> &tvec, const Approx &cr) {
  if (verifyVectorApproxCompatibility(tvec, cr) == false) {
    return false;
  }
  const size_t npts = tvec.size();
  const std::vector<double> crval = cr.getValues();
  const double cr_margin = cr.getMargin();
  for (size_t i = 0; i < npts; i++) {
    if (tvec[i] <= crval[i] + cr_margin) {
      return false;
    }
  }
  return true;
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool operator>(const Approx &cr, const std::vector<T> &tvec) {
  if (verifyVectorApproxCompatibility(tvec, cr) == false) {
    return false;
  }
  const size_t npts = tvec.size();
  const std::vector<double> crval = cr.getValues();
  const double cr_margin = cr.getMargin();
  for (size_t i = 0; i < npts; i++) {
    if (crval[i] - cr_margin <= tvec[i]) {
      return false;
    }
  }
  return true;
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool operator<(const std::vector<T> &tvec, const Approx &cr) {
  if (verifyVectorApproxCompatibility(tvec, cr) == false) {
    return false;
  }
  const size_t npts = tvec.size();
  const std::vector<double> crval = cr.getValues();
  const double cr_margin = cr.getMargin();
  for (size_t i = 0; i < npts; i++) {
    if (tvec[i] >= crval[i] - cr_margin) {
      return false;
    }
  }
  return true;
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool operator<(const Approx &cr, const std::vector<T> &tvec) {
  if (verifyVectorApproxCompatibility(tvec, cr) == false) {
    return false;
  }
  const size_t npts = tvec.size();
  const std::vector<double> crval = cr.getValues();
  const double cr_margin = cr.getMargin();
  for (size_t i = 0; i < npts; i++) {
    if (crval[i] + cr_margin >= tvec[i]) {
      return false;
    }
  }
  return true;
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool operator>=(const std::vector<T> &tvec, const Approx &cr) {
  if (verifyVectorApproxCompatibility(tvec, cr) == false) {
    return false;
  }
  const size_t npts = tvec.size();
  const std::vector<double> crval = cr.getValues();
  const double cr_margin = cr.getMargin();
  for (size_t i = 0; i < npts; i++) {
    if (tvec[i] < crval[i] - cr_margin) {
      return false;
    }
  }
  return true;
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool operator>=(const Approx &cr, const std::vector<T> &tvec) {
  if (verifyVectorApproxCompatibility(tvec, cr) == false) {
    return false;
  }
  const size_t npts = tvec.size();
  const std::vector<double> crval = cr.getValues();
  const double cr_margin = cr.getMargin();
  for (size_t i = 0; i < npts; i++) {
    if (crval[i] + cr_margin < tvec[i]) {
      return false;
    }
  }
  return true;
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool operator<=(const std::vector<T> &tvec, const Approx &cr) {
  if (verifyVectorApproxCompatibility(tvec, cr) == false) {
    return false;
  }
  const size_t npts = tvec.size();
  const std::vector<double> crval = cr.getValues();
  const double cr_margin = cr.getMargin();
  for (size_t i = 0; i < npts; i++) {
    if (tvec[i] > crval[i] + cr_margin) {
      return false;
    }
  }
  return true;
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool operator<=(const Approx &cr, const std::vector<T> &tvec) {
  if (verifyVectorApproxCompatibility(tvec, cr) == false) {
    return false;
  }
  const size_t npts = tvec.size();
  const std::vector<double> crval = cr.getValues();
  const double cr_margin = cr.getMargin();
  for (size_t i = 0; i < npts; i++) {
    if (crval[i] - cr_margin > tvec[i]) {
      return false;
    }
  }
  return true;
}

} // namespace testing
} // namespace stormm
