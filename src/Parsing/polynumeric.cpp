#include <cstring>
#include "copyright.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "polynumeric.h"

namespace stormm {
namespace parse {

using data_types::getStormmScalarTypeName;
using data_types::getStormmHpcVectorTypeName;

//-------------------------------------------------------------------------------------------------
std::vector<double> doubleFromPolyNumeric(const std::vector<PolyNumeric> &values) {
  const size_t vsize = values.size();
  std::vector<double> result(vsize, 0.0);
  for (size_t i = 0; i < vsize; i++) {
    result[i] = values[i].d;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> intFromPolyNumeric(const std::vector<PolyNumeric> &values) {
  const size_t vsize = values.size();
  std::vector<int> result(vsize, 0);
  for (size_t i = 0; i < vsize; i++) {
    result[i] = values[i].i;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<llint> llintFromPolyNumeric(const std::vector<PolyNumeric> &values) {
  const size_t vsize = values.size();
  std::vector<llint> result(vsize, 0LL);
  for (size_t i = 0; i < vsize; i++) {
    result[i] = values[i].lli;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> uintFromPolyNumeric(const std::vector<PolyNumeric> &values) {
  const size_t vsize = values.size();
  std::vector<uint> result(vsize, 0);
  for (size_t i = 0; i < vsize; i++) {
    result[i] = values[i].ui;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<ullint> ullintFromPolyNumeric(const std::vector<PolyNumeric> &values) {
  const size_t vsize = values.size();
  std::vector<ullint> result(vsize, 0LLU);
  for (size_t i = 0; i < vsize; i++) {
    result[i] = values[i].ulli;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<char4> char4FromPolyNumeric(const std::vector<PolyNumeric> &values) {
  const size_t vsize = values.size();
  char4 tmp = {'\0', '\0', '\0', '\0'};
  std::vector<char4> result(vsize, tmp);
  for (size_t i = 0; i < vsize; i++) {
    result[i] = values[i].c4;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<PolyNumeric> polyNumericVector(const std::vector<char4> &values) {
  std::vector<PolyNumeric> result;
  const size_t vsize = values.size();
  result.resize(vsize);
  for (size_t i = 0; i < vsize; i++) {
    result[i].c4.x = values[i].x;
    result[i].c4.y = values[i].y;
    result[i].c4.z = values[i].z;
    result[i].c4.w = values[i].w;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string nameNumericalType(NumberFormat cform) {
  switch (cform) {
  case NumberFormat::SCIENTIFIC:
  case NumberFormat::STANDARD_REAL:
    return getStormmScalarTypeName<double>();
  case NumberFormat::INTEGER:
    return getStormmScalarTypeName<int>();
  case NumberFormat::LONG_LONG_INTEGER:
    return getStormmScalarTypeName<llint>();
  case NumberFormat::UNSIGNED_INTEGER:
    return getStormmScalarTypeName<ulint>();
  case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
    return getStormmScalarTypeName<ullint>();
  case NumberFormat::CHAR4:
    return getStormmHpcVectorTypeName<char4>();
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string nameNumberFormat(NumberFormat cform) {
  switch (cform) {
  case NumberFormat::SCIENTIFIC:
    return std::string("SCIENTIFIC");
  case NumberFormat::STANDARD_REAL:
    return std::string("STANDARD_REAL");
  case NumberFormat::INTEGER:
    return std::string("INTEGER");
  case NumberFormat::LONG_LONG_INTEGER:
    return std::string("LONG_LONG_INTEGER");
  case NumberFormat::UNSIGNED_INTEGER:
    return std::string("UNSIGNED_INTEGER");
  case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
    return std::string("UNSIGNED_LONG_LONG_INTEGER");
  case NumberFormat::CHAR4:
    return std::string("CHAR4");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
PolyNumeric extractFormattedNumber(const char* a, const NumberFormat cform, const int read_begin,
                                   const int len) {
  PolyNumeric pn;

  // Pre-allocate a buffer and fill it with the text of interest
  const int width = (len > 0) ? len : strlen(a) - read_begin;
  std::vector<char> buffer(width + 1, '\0');
  for (int i = 0; i < width; i++) {
    buffer[i] = a[read_begin + i];
  }
  switch (cform) {
  case NumberFormat::SCIENTIFIC:
  case NumberFormat::STANDARD_REAL:
    pn.d = strtod(buffer.data(), nullptr);
    break;
  case NumberFormat::INTEGER:
    pn.i = strtol(buffer.data(), nullptr, 10);
    break;
  case NumberFormat::LONG_LONG_INTEGER:
    pn.lli = strtoll(buffer.data(), nullptr, 10);
    break;
  case NumberFormat::UNSIGNED_INTEGER:
    pn.ui = strtoul(buffer.data(), nullptr, 10);
    break;
  case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
    pn.ulli = strtoull(buffer.data(), nullptr, 10);
    break;
  case NumberFormat::CHAR4:
    {
      int i = 0;
      while (i < width && buffer[i] == ' ') {
        i++;
      }
      if (i == 0) {
        pn.c4.x = buffer[0];
        pn.c4.y = buffer[1];
        pn.c4.z = buffer[2];
        pn.c4.w = buffer[3];
      }
      else {
        pn.c4.x = (i < 4) ? buffer[i    ] : ' ';
        pn.c4.y = (i < 3) ? buffer[i + 1] : ' ';
        pn.c4.z = (i < 2) ? buffer[i + 2] : ' ';
        pn.c4.w = (i < 1) ? buffer[i + 3] : ' ';
      }
    }
    break;
  }
  return pn;
}

} // namespace parse
} // namespace stormm
