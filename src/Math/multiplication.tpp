// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
template <typename T> double logProduct(const T* values, const size_t length) {
  double esum = 0.0;
  for (size_t i = 0LLU; i < length; i++) {
    if (values[i] <= (T)(0)) {
      rtErr("The logarithm of a negative number, or zero, is undefined.  " +
            realToString(values[i], 11, 4, NumberFormat::STANDARD_REAL) + " was encountered in "
            "position " + std::to_string(i) + ".", "logProduct");
    }
    esum += log(values[i]);
  }
  return esum;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double logProduct(const std::vector<T> &values) {
  return logProduct(values.data(), values.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double logProduct(const Hybrid<T> &values) {
  return logProduct(values.data(), values.size());
}

//-------------------------------------------------------------------------------------------------
template <typename Tprod, typename Tbase> Tprod seriesProduct(const Tbase* va,
                                                              const size_t length) {
  if (isScalarType<Tbase>() == false) {
    rtErr("Data type " + getStormmScalarTypeName<Tbase>() + " is not suitable for computing the "
          "product series of scalar values.", "productSeries");
  }
  if (isScalarType<Tprod>() == false) {
    rtErr("Data type " + getStormmScalarTypeName<Tprod>() + " is not suitable for representing "
          "the product of a series of scalar numbers.", "productSeries");
  }
  if (isSignedIntegralScalarType<Tbase>() || isUnsignedIntegralScalarType<Tbase>()) {
    const double base_e_limit = (sizeof(Tprod) == 4) ? log(2.0) * 31 : log(2.0) * 63;
    double log_p = 0.0;
    for (size_t i = 0LLU; i < length; i++) {
      if (va[i] == 0.0) {
        return static_cast<Tprod>(0);
      }
      else if (va[i] < 0.0) {
        log_p += log(-va[i]);
      }
      else {
        log_p += log(va[i]);
      }
    }
    if (log_p >= base_e_limit) {
      rtErr("The product of a series of " + std::to_string(length) + " integers comes out greater "
            "than the maximum representable value of the chosen format (" +
            getStormmScalarTypeName<Tprod>() + ").", "seriesProduct");
    }
    Tprod result = 1;
    for (size_t i = 0LLU; i < length; i++) {
      result *= va[i];
    }
    return result;
  }
  else if (isFloatingPointScalarType<Tbase>()) {
    if (isFloatingPointScalarType<Tprod>() == false) {
      rtErr("Computing the product of a series of floating point numbers (" +
            getStormmScalarTypeName<Tbase>() + ") as an integral type (" +
            getStormmScalarTypeName<Tprod>() + ") is unsafe and not permitted.", "seriesProduct");
    }

    // Compute the logarithm of the product to ensure that it stays within the limits of the
    // hardware.
    double log_p = 0.0;
    const double base_e_limit = (sizeof(Tprod) == 4) ? log(2.0) * 128.0 : log(2.0) * 1024.0;
    bool violation = false;
    double sign_mult = 1.0;
    for (size_t i = 0LLU; i < length; i++) {
      if (va[i] == 0.0) {
        return 0.0;
      }
      else if (va[i] < 0.0) {
        log_p += log(-va[i]);
        sign_mult *= -1.0;
      }
      else {
        log_p += log(va[i]);
      }
      violation = (violation || (fabs(log_p) >= base_e_limit));
    }
    if (violation) {
      if (log_p < base_e_limit) {
        return sign_mult * exp(log_p);
      }
      else {
        rtErr("The product of a series of " + std::to_string(length) + " real-valued numbers "
              "comes out greater than the maximum representable format in " +
              getStormmScalarTypeName<Tprod>() + " numbers.  Choose a longer format, or use the "
              "logProduct() function to calculate the natural logarithm of the product.",
              "seriesProduct");
      }
    }
    else {
      Tprod result = 1.0;
      for (size_t i = 0LLU; i < length; i++) {
        result *= va[i];
      }
      return result;
    }
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tprod, typename Tbase> Tprod seriesProduct(const std::vector<Tbase> &va) {
  return seriesProduct<Tprod, Tbase>(va.data(), va.size());
}

//-------------------------------------------------------------------------------------------------
template <typename Tprod, typename Tbase> Tprod seriesProduct(const Hybrid<Tbase> &va) {
  return seriesProduct<Tprod, Tbase>(va.data(), va.size());
}

} // namespace stmath
} // namespace stormm
