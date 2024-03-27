#include <cmath>
#include "copyright.h"
#include "approx.h"

namespace stormm {
namespace testing {

//-------------------------------------------------------------------------------------------------
Approx::Approx(const double value_in, ComparisonType style_in, const double tol_in) :
    Approx(std::vector<double>(1, value_in), style_in, tol_in)
{}

//-------------------------------------------------------------------------------------------------
Approx::Approx(const double value_in, const double tol_in, ComparisonType style_in) :
    Approx(std::vector<double>(1, value_in), style_in, tol_in)
{}

//-------------------------------------------------------------------------------------------------
int Approx::size() const {
  return values.size();
}

//-------------------------------------------------------------------------------------------------
double Approx::getValue() const {
  if (values.size() == 1) {
    return values[0];
  }
  else {
    rtErr("A single value was requested for an approximate comparison containing a vector of " +
          std::to_string(values.size()) + " values.  Use getValues() to retrieve them all.");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<double> Approx::getValues() const {
  return values;
}

//-------------------------------------------------------------------------------------------------
ComparisonType Approx::getStyle() const {
  return style;
}

//-------------------------------------------------------------------------------------------------
double Approx::getMargin() const {
  return dtol;
}

//-------------------------------------------------------------------------------------------------
double Approx::getTolerance() const {
  return dtol;
}

//-------------------------------------------------------------------------------------------------
double Approx::getTol() const {
  return dtol;
}

//-------------------------------------------------------------------------------------------------
void Approx::setValue(const double value_in) {
  values.resize(1);
  values[0] = value_in;
  values.shrink_to_fit();
}

//-------------------------------------------------------------------------------------------------
void Approx::setValues(const std::vector<double> &values_in) {
  values = values_in;
}

//-------------------------------------------------------------------------------------------------
void Approx::setValue(const double value_in, const size_t index) {
  if (index >= values.size()) {
    rtErr("Index " + std::to_string(index) + " is invalid for an approximate comparison with " +
          std::to_string(values.size()) + " entries.", "Approx", "setValue");
  }
  values[index] = value_in;
}

//-------------------------------------------------------------------------------------------------
void Approx::setMargin(const double dtol_in) {
  dtol = dtol_in;
}

//-------------------------------------------------------------------------------------------------
void Approx::setTolerance(const double dtol_in) {
  dtol = dtol_in;
}

//-------------------------------------------------------------------------------------------------
void Approx::setTol(const double dtol_in) {
  dtol = dtol_in;
}

//-------------------------------------------------------------------------------------------------
Approx Approx::margin(const double dtol_in) const {
  return Approx(values, style, dtol_in);
}

//-------------------------------------------------------------------------------------------------
Approx Approx::tolerance(const double dtol_in) const {
  return Approx(values, style, dtol_in);
}

//-------------------------------------------------------------------------------------------------
Approx Approx::tol(const double dtol_in) const {
  return Approx(values, style, dtol_in);
}

//-------------------------------------------------------------------------------------------------
bool Approx::test(const double test_value) const {
  if (values.size() != 1) {
    return false;
  }
  switch (style) {
  case ComparisonType::ABSOLUTE:
  case ComparisonType::MEAN_UNSIGNED_ERROR:
    return (std::abs(test_value - values[0]) <= dtol);
  case ComparisonType::RELATIVE:
  case ComparisonType::RELATIVE_RMS_ERROR:
    if (std::abs(values[0]) > constants::tiny) {
      return (std::abs((test_value - values[0]) / values[0]) <= dtol);
    }
    else {
      return (std::abs((test_value - values[0]) / constants::tiny) <= dtol);
    }
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool operator==(const double d, const Approx &cr) {
  return cr.test(d);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const Approx &cr, const double d) {
  return cr.test(d);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const double d, const Approx &cr) {
  return (cr.test(d) == false);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const Approx &cr, const double d) {
  return (cr.test(d) == false);
}

//-------------------------------------------------------------------------------------------------
bool operator>(const double d, const Approx &cr) {
  return (d > cr.getValue() + cr.getMargin());
}

//-------------------------------------------------------------------------------------------------
bool operator>(const Approx &cr, const double d) {
  return (cr.getValue() - cr.getMargin() > d);
}

//-------------------------------------------------------------------------------------------------
bool operator<(const double d, const Approx &cr) {
  return (d < cr.getValue() - cr.getMargin());
}

//-------------------------------------------------------------------------------------------------
bool operator<(const Approx &cr, const double d) {
  return (cr.getValue() + cr.getMargin() < d);
}

//-------------------------------------------------------------------------------------------------
bool operator>=(const double d, const Approx &cr) {
  return (d >= cr.getValue() - cr.getMargin());
}

//-------------------------------------------------------------------------------------------------
bool operator>=(const Approx &cr, const double d) {
  return (cr.getValue() + cr.getMargin() >= d);
}

//-------------------------------------------------------------------------------------------------
bool operator<=(const double d, const Approx &cr) {
  return (d <= cr.getValue() + cr.getMargin());
}

//-------------------------------------------------------------------------------------------------
bool operator<=(const Approx &cr, const double d) {
  return (cr.getValue() - cr.getMargin() <= d);
}

} // namespace testing
} // namespace stormm
