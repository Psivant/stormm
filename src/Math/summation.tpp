// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace stmath {

using data_types::isScalarType;
using data_types::isFloatingPointScalarType;
using data_types::getStormmScalarTypeName;

//-------------------------------------------------------------------------------------------------
template <typename TBase>
void prefixSumInPlace(TBase* vdata, const size_t n_elements, const PrefixSumType style,
                      const char* caller) {
  TBase sum = 0;
  if (isFloatingPointScalarType<TBase>()) {
    double lfsum = 0.0;
    switch (style) {
    case PrefixSumType::EXCLUSIVE:
      for (size_t i = 0; i < n_elements; i++) {
        const TBase tmp_sum = vdata[i];
        vdata[i] = sum;
        sum += tmp_sum;
        lfsum += static_cast<double>(tmp_sum);
      }
      break;
    case PrefixSumType::INCLUSIVE:
      for (size_t i = 0; i < n_elements; i++) {
        sum += vdata[i];
        lfsum += vdata[i];
        vdata[i] = sum;
      }
      break;
    }
    if (lfsum != Approx(sum, ComparisonType::RELATIVE, constants::small)) {
      const std::string tbase_name = isScalarType<TBase>() ?
                                     getStormmScalarTypeName<TBase>() :
                                     std::string(std::type_index(typeid(TBase)).name());
      const std::string callfunc = (caller == nullptr) ? "" :
                                                         ".  Called by " + std::string(caller);
      rtErr("Overflow of numerical format.  Summation of a " + std::to_string(n_elements) +
            "-element series as double produces " + std::to_string(lfsum) +
            ", whereas the array's base type " + tbase_name + " records " + std::to_string(sum) +
            callfunc + ".", "prefixSumInPlace");
    }
  }
  else {
    llint llsum = 0;
    switch (style) {
    case PrefixSumType::EXCLUSIVE:
      for (size_t i = 0; i < n_elements; i++) {
        const TBase tmp_sum = vdata[i];
        vdata[i] = sum;
        sum += tmp_sum;
        llsum += static_cast<llint>(tmp_sum);
      }
      break;
    case PrefixSumType::INCLUSIVE:
      for (size_t i = 0; i < n_elements; i++) {
        sum += vdata[i];
        llsum += vdata[i];
        vdata[i] = sum;
      }
      break;
    }
    if (llsum != Approx(sum, ComparisonType::RELATIVE, constants::small)) {
      const std::string tbase_name = isScalarType<TBase>() ?
                                     getStormmScalarTypeName<TBase>() :
                                     std::string(std::type_index(typeid(TBase)).name());
      const std::string callfunc = (caller == nullptr) ? "" :
                                                         ".  Called by " + std::string(caller);
      rtErr("Overflow of numerical format.  Summation of a " + std::to_string(n_elements) +
            "-element series as llint produces " + std::to_string(llsum) +
            ", whereas the array's base type " + tbase_name + " records " + std::to_string(sum) +
            callfunc + ".", "prefixSumInPlace");
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename TBase>
void prefixSumInPlace(std::vector<TBase> *v, const PrefixSumType style, const char* caller) {
  prefixSumInPlace(v->data(), v->size(), style, caller);
}

//-------------------------------------------------------------------------------------------------
template <typename TBase>
void prefixSumInPlace(Hybrid<TBase> *v, const PrefixSumType style, const char* caller) {
  prefixSumInPlace(v->data(), v->size(), style, caller);
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sum(const TBase* v, const size_t vlen) {
  TSum total = static_cast<TSum>(0);
  for (size_t i = 0; i < vlen; i++) {
    total += static_cast<TSum>(v[i]);
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sum(const std::vector<TBase> &v) {
  const size_t vlen = v.size();
  TSum total = static_cast<TSum>(0);
  for (size_t i = 0; i < vlen; i++) {
    total += static_cast<TSum>(v[i]);
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase> TSum sum(const Hybrid<TBase> &hb) {
  const size_t hblen = hb.size();
  TSum total = static_cast<TSum>(0);
  const TBase* hbptr = hb.data();
  for (size_t i = 0; i < hblen; i++) {
    total += static_cast<TSum>(hbptr[i]);
  }
  return total;
}
  
//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sumTuple2(const TBase* v, const size_t vlen) {
  TSum total = { 0, 0 };
  for (size_t i = 0; i < vlen; i++) {
    total.x += v[i].x;
    total.y += v[i].y;
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sumTuple2(const std::vector<TBase> &v) {
  const size_t vlen = v.size();  
  TSum total = { 0, 0 };
  for (size_t i = 0; i < vlen; i++) {
    total.x += v[i].x;
    total.y += v[i].y;
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sumTuple2(const Hybrid<TBase> &hb) {
  const size_t hblen = hb.size();
  TSum total = { 0, 0 };
  const TBase* hbptr = hb.data();
  for (size_t i = 0; i < hblen; i++) {
    total.x += static_cast<TSum>(hbptr[i].x);
    total.y += static_cast<TSum>(hbptr[i].y);
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sumTuple3(const TBase* v, const size_t vlen) {
  TSum total = { 0, 0, 0 };
  for (size_t i = 0; i < vlen; i++) {
    total.x += v[i].x;
    total.y += v[i].y;
    total.z += v[i].z;
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sumTuple3(const std::vector<TBase> &v) {
  const size_t vlen = v.size();  
  TSum total = { 0, 0, 0 };
  for (size_t i = 0; i < vlen; i++) {
    total.x += v[i].x;
    total.y += v[i].y;
    total.z += v[i].z;
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sumTuple3(const Hybrid<TBase> &hb) {
  const size_t hblen = hb.size();
  TSum total = { 0, 0, 0 };
  const TBase* hbptr = hb.data();
  for (size_t i = 0; i < hblen; i++) {
    total.x += static_cast<TSum>(hbptr[i].x);
    total.y += static_cast<TSum>(hbptr[i].y);
    total.z += static_cast<TSum>(hbptr[i].z);
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sumTuple4(const TBase* v, const size_t vlen) {
  TSum total = { 0, 0, 0, 0 };
  for (size_t i = 0; i < vlen; i++) {
    total.x += v[i].x;
    total.y += v[i].y;
    total.z += v[i].z;
    total.w += v[i].w;
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sumTuple4(const std::vector<TBase> &v) {
  const size_t vlen = v.size();  
  TSum total = { 0, 0, 0, 0 };
  for (size_t i = 0; i < vlen; i++) {
    total.x += v[i].x;
    total.y += v[i].y;
    total.z += v[i].z;
    total.w += v[i].w;
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename TSum, typename TBase>
TSum sumTuple4(const Hybrid<TBase> &hb) {
  const size_t hblen = hb.size();
  TSum total = { 0, 0, 0, 0 };
  const TBase* hbptr = hb.data();
  for (size_t i = 0; i < hblen; i++) {
    total.x += static_cast<TSum>(hbptr[i].x);
    total.y += static_cast<TSum>(hbptr[i].y);
    total.z += static_cast<TSum>(hbptr[i].z);
    total.w += static_cast<TSum>(hbptr[i].w);
  }
  return total;
}

//-------------------------------------------------------------------------------------------------
template <typename Ta, typename Tb>
std::vector<Ta> sum(const std::vector<Ta> &va, const std::vector<Tb> &vb, const double afac,
                    const double bfac) {

  // Check that all parameters have the same length
  const size_t vlen = va.size();
  if (vlen != vb.size()) {
    rtErr("Vectors of differing sizes (" + std::to_string(vlen) + ", " +
          std::to_string(vb.size()) + ") cannot be summed in element-wise fashion.", "sum");
  }
  std::vector<Ta> result(vlen);
  if (afac == 1.0 && bfac == 1.0) {
    for (size_t i = 0; i < vlen; i++) {
      result[i] = va[i] + static_cast<Ta>(vb[i]);
    }
  }
  else {
    for (size_t i = 0; i < vlen; i++) {
      result[i] = static_cast<Ta>((static_cast<double>(va[i]) * afac) +
                                  (static_cast<double>(vb[i]) * bfac));
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Ta, typename Tb, typename Tc>
std::vector<Ta> sum(const std::vector<Ta> &va, const std::vector<Tb> &vb,
                    const std::vector<Tc> &vc, const double afac, const double bfac,
                    const double cfac) {

  // Check that all parameters have the same length
  const size_t vlen = va.size();
  if (vlen != vb.size() || vlen != vc.size()) {
    rtErr("Vectors of differing sizes (" + std::to_string(vlen) + ", " +
          std::to_string(vb.size()) + ", " + std::to_string(vc.size()) + ") cannot "
          "be summed in element-wise fashion.", "sum");
  }
  std::vector<Ta> result(vlen);
  if (afac == 1.0 && bfac == 1.0 && cfac == 1.0) {
    for (size_t i = 0; i < vlen; i++) {
      result[i] = va[i] + static_cast<Ta>(vb[i]) + static_cast<Ta>(vc[i]);
    }
  }
  else {
    for (size_t i = 0; i < vlen; i++) {
      result[i] = static_cast<Ta>((static_cast<double>(va[i]) * afac) +
                                  (static_cast<double>(vb[i]) * bfac) +
                                  (static_cast<double>(vc[i]) * cfac));
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Ta, typename Tb, typename Tc, typename Td>
std::vector<Ta> sum(const std::vector<Ta> &va, const std::vector<Tb> &vb,
                    const std::vector<Tc> &vc, const std::vector<Td> &vd, const double afac,
                    const double bfac, const double cfac, const double dfac) {

  // Check that all parameters have the same length
  const size_t vlen = va.size();
  if (vlen != vb.size() || vlen != vc.size() || vlen != vd.size()) {
    rtErr("Vectors of differing sizes (" + std::to_string(vlen) + ", " +
          std::to_string(vb.size()) + ", " + std::to_string(vc.size()) + ", " +
          std::to_string(vd.size()) + ") cannot be summed in element-wise fashion.", "sum");
  }
  std::vector<Ta> result(vlen);
  if (afac == 1.0 && bfac == 1.0 && cfac == 1.0 && dfac == 1.0) {
    for (size_t i = 0; i < vlen; i++) {
      result[i] = va[i] + static_cast<Ta>(vb[i]) + static_cast<Ta>(vc[i]) + static_cast<Ta>(vd[i]);
    }
  }
  else {
    for (size_t i = 0; i < vlen; i++) {
      result[i] = static_cast<Ta>((static_cast<double>(va[i]) * afac) +
                                  (static_cast<double>(vb[i]) * bfac) +
                                  (static_cast<double>(vc[i]) * cfac) +
                                  (static_cast<double>(vd[i]) * dfac));
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Ta, typename Tb>
void sum(Ta* va, const Tb* vb, const size_t length, const double afac, const double bfac) {
  if (afac == 1.0 && bfac == 1.0) {
    for (size_t i = 0; i < length; i++) {
      va[i] += static_cast<Ta>(vb[i]);
    }
  }
  else {
    for (size_t i = 0; i < length; i++) {
      va[i] = static_cast<Ta>((static_cast<double>(va[i]) * afac) +
                              (static_cast<double>(vb[i]) * bfac));
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Ta, typename Tb, typename Tc>
void sum(Ta* va, const Tb* vb, const Tc* vc, const size_t length, const double afac,
         const double bfac, const double cfac) {
  if (afac == 1.0 && bfac == 1.0 && cfac == 1.0) {
    for (size_t i = 0; i < length; i++) {
      va[i] += static_cast<Ta>(vb[i]) + static_cast<Ta>(vc[i]);
    }
  }
  else {
    for (size_t i = 0; i < length; i++) {
      va[i] = static_cast<Ta>((static_cast<double>(va[i]) * afac) +
                              (static_cast<double>(vb[i]) * bfac) +
                              (static_cast<double>(vc[i]) * cfac));
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Ta, typename Tb, typename Tc, typename Td>
void sum(Ta* va, const Tb* vb, const Tc* vc, const Td* vd, const size_t length, const double afac,
         const double bfac, const double cfac, const double dfac) {
  if (afac == 1.0 && bfac == 1.0 && cfac == 1.0 && dfac == 1.0) {
    for (size_t i = 0; i < length; i++) {
      va[i] += static_cast<Ta>(vb[i]) + static_cast<Ta>(vc[i]) + static_cast<Ta>(vd[i]);
    }
  }
  else {
    for (size_t i = 0; i < length; i++) {
      va[i] = static_cast<Ta>((static_cast<double>(va[i]) * afac) +
                              (static_cast<double>(vb[i]) * bfac) +
                              (static_cast<double>(vc[i]) * cfac) +
                              (static_cast<double>(vd[i]) * dfac));
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Ta, typename Tb>
void sum(std::vector<Ta> *va, const std::vector<Tb> &vb, const double afac, const double bfac) {

  // Check that all parameters have the same length
  const size_t vlen = va->size();
  if (vlen != vb.size()) {
    rtErr("Vectors of differing sizes (" + std::to_string(vlen) + ", " +
          std::to_string(vb.size()) + ") cannot be summed in element-wise fashion.", "sum");
  }
  sum(va->data(), vb.data(), vlen, afac, bfac);
}

//-------------------------------------------------------------------------------------------------
template <typename Ta, typename Tb, typename Tc>
void sum(std::vector<Ta> *va, const std::vector<Tb> &vb, const std::vector<Tc> &vc,
         const double afac, const double bfac, const double cfac) {

  // Check that all parameters have the same length
  const size_t vlen = va->size();
  if (vlen != vb.size() || vlen != vc.size()) {
    rtErr("Vectors of differing sizes (" + std::to_string(vlen) + ", " +
          std::to_string(vb.size()) + ", " + std::to_string(vc.size()) + ") cannot be summed in "
          "element-wise fashion.", "sum");
  }
  sum(va->data(), vb.data(), vc.data(), vlen, afac, bfac, cfac);
}

//-------------------------------------------------------------------------------------------------
template <typename Ta, typename Tb, typename Tc, typename Td>
void sum(std::vector<Ta> *va, const std::vector<Tb> &vb, const std::vector<Tc> &vc,
         const std::vector<Td> &vd, const double afac, const double bfac, const double cfac,
         const double dfac) {

  // Check that all parameters have the same length
  const size_t vlen = va->size();
  if (vlen != vb.size() || vlen != vc.size() || vlen != vd.size()) {
    rtErr("Vectors of differing sizes (" + std::to_string(vlen) + ", " +
          std::to_string(vb.size()) + ", " + std::to_string(vc.size()) + ", " +
          std::to_string(vd.size()) + ") cannot be summed in element-wise fashion.", "sum");
  }
  sum(va->data(), vb.data(), vc.data(), vd.data(), vlen, afac, bfac, cfac, dfac);
}

//-------------------------------------------------------------------------------------------------
template <typename Ta, typename Tb>
void sum(Hybrid<Ta> *va, const Hybrid<Tb> &vb, const double afac, const double bfac) {

  // Check that all parameters have the same length
  const size_t vlen = va->size();
  if (vlen != vb.size()) {
    rtErr("Vectors of differing sizes (" + std::to_string(vlen) + ", " +
          std::to_string(vb.size()) + ") cannot be summed in element-wise fashion.", "sum");
  }
  sum(va->data(), vb.data(), vlen, afac, bfac);
}

//-------------------------------------------------------------------------------------------------
template <typename Ta, typename Tb, typename Tc>
void sum(Hybrid<Ta> *va, const Hybrid<Tb> &vb, const Hybrid<Tc> &vc, const double afac,
         const double bfac, const double cfac) {

  // Check that all parameters have the same length
  const size_t vlen = va->size();
  if (vlen != vb.size() || vlen != vc.size()) {
    rtErr("Vectors of differing sizes (" + std::to_string(vlen) + ", " +
          std::to_string(vb.size()) + ", " + std::to_string(vc.size()) + ") cannot be summed in "
          "element-wise fashion.", "sum");
  }
  sum(va->data(), vb.data(), vc.data(), vlen, afac, bfac, cfac);
}

//-------------------------------------------------------------------------------------------------
template <typename Ta, typename Tb, typename Tc, typename Td>
void sum(Hybrid<Ta> *va, const Hybrid<Tb> &vb, const Hybrid<Tc> &vc, const Hybrid<Td> &vd,
         const double afac, const double bfac, const double cfac, const double dfac) {

  // Check that all parameters have the same length
  const size_t vlen = va->size();
  if (vlen != vb.size() || vlen != vc.size() || vlen != vd.size()) {
    rtErr("Vectors of differing sizes (" + std::to_string(vlen) + ", " +
          std::to_string(vb.size()) + ", " + std::to_string(vc.size()) + ", " +
          std::to_string(vd.size()) + ") cannot be summed in element-wise fashion.", "sum");
  }
  sum(va->data(), vb.data(), vc.data(), vd.data(), vlen, afac, bfac, cfac, dfac);
}

} // namespace stmath
} // namespace stormm
