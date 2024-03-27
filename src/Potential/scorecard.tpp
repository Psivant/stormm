// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename T> T ScoreCard::getEnergyScalingFactor() const {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  T result;
  if (ct == float_type_index) {
    result = nrg_scale_f;
  }
  else if (ct == double_type_index) {
    result = nrg_scale_lf;
  }
  else {
    rtErr("No scaling factor is defined for data in format " +
          std::string(std::type_index(typeid(T)).name()) + ".", "ScoreCard",
          "getEnergyScalingFactor");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> T ScoreCard::getInverseEnergyScalingFactor() const {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  T result;
  if (ct == float_type_index) {
    result = inverse_nrg_scale_f;
  }
  else if (ct == double_type_index) {
    result = inverse_nrg_scale_lf;
  }
  else {
    rtErr("No scaling factor is defined for data in format " +
          std::string(std::type_index(typeid(T)).name()) + ".", "ScoreCard",
          "getInverseEnergyScalingFactor");
  }
  return result;
}

} // namespace energy
} // namespace stormm
