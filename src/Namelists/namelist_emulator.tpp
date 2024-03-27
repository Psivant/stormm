// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace namelist {

//-------------------------------------------------------------------------------------------------
template <typename T>
void NamelistEmulator::assignVariable(T *var, const double mult, const std::string &keyword_query,
                                      const int index) const {
  if (getKeywordStatus(keyword_query) != InputStatus::MISSING) {
    if (isSignedIntegralScalarType<T>() || isUnsignedIntegralScalarType<T>()) {
      *var = getIntValue(keyword_query, index) * static_cast<T>(mult);
    }
    else if (isFloatingPointScalarType<T>()) {
      *var = getRealValue(keyword_query, index) * mult;
    }
    else {
      rtErr("No conversion is available for data type " +
            std::string(std::type_index(typeid(T)).name()) + ".", "NamelistEmulator",
            "assignVariable");
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void NamelistEmulator::assignVariable(T *var, const std::string &keyword_query,
                                      const int index) const {
  assignVariable<T>(var, 1.0, keyword_query, index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void NamelistEmulator::assignVariable(T *var_x, T *var_y, T *var_z, const double mult,
                                      const std::string &keyword_query, const int index) const {
  if (getKeywordStatus(keyword_query) != InputStatus::MISSING) {
    if (isSignedIntegralScalarType<T>() || isUnsignedIntegralScalarType<T>()) {
      *var_x = getIntValue(keyword_query, index) * static_cast<T>(mult);
      *var_y = *var_x;
      *var_z = *var_x;
    }
    else if (isFloatingPointScalarType<T>()) {
      *var_x = getRealValue(keyword_query, index) * mult;
      *var_y = *var_x;
      *var_z = *var_x;
    }
    else {
      rtErr("Triplicate conversion is not available for data type " +
            std::string(std::type_index(typeid(T)).name()) + ".", "NamelistEmulator",
            "assignVariable");
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void NamelistEmulator::assignVariable(T *var_x, T *var_y, T *var_z,
                                      const std::string &keyword_query, const int index) const {
  assignVariable<T>(var_x, var_y, var_z, 1.0, keyword_query, index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void NamelistEmulator::assignVariable(T *var, const double mult, const std::string &keyword_query,
                                      const std::string &sub_key, const int index) const {
  if (getKeywordStatus(keyword_query, sub_key, index) != InputStatus::MISSING) {
    if (isSignedIntegralScalarType<T>() || isUnsignedIntegralScalarType<T>()) {
      *var = getIntValue(keyword_query, sub_key, index) * static_cast<T>(mult);
    }
    else if (isFloatingPointScalarType<T>()) {
      *var = getRealValue(keyword_query, sub_key, index) * mult;
    }
    else {
      rtErr("No conversion is available for data type " +
            std::string(std::type_index(typeid(T)).name()) + ".", "NamelistEmulator",
            "assignVariable");
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void NamelistEmulator::assignVariable(T *var, const std::string &keyword_query,
                                      const std::string &sub_key, const int index) const {
  assignVariable<T>(var, 1.0, keyword_query, sub_key, index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void NamelistEmulator::assignVariable(T *var_x, T *var_y, T *var_z, const double mult,
                                      const std::string &keyword_query, const std::string &sub_key,
                                      const int index) const {
  if (getKeywordStatus(keyword_query, sub_key, index) != InputStatus::MISSING) {
    if (isSignedIntegralScalarType<T>() || isUnsignedIntegralScalarType<T>()) {
      *var_x = getIntValue(keyword_query, sub_key, index) * static_cast<T>(mult);
      *var_y = *var_x;
      *var_z = *var_x;
    }
    else if (isFloatingPointScalarType<T>()) {
      *var_x = getRealValue(keyword_query, sub_key, index) * mult;
      *var_y = *var_x;
      *var_z = *var_x;
    }
    else {
      rtErr("Triplicate conversion is not available for data type " +
            std::string(std::type_index(typeid(T)).name()) + ".", "NamelistEmulator",
            "assignVariable");
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void NamelistEmulator::assignVariable(T *var_x, T *var_y, T *var_z,
                                      const std::string &keyword_query, const std::string &sub_key,
                                      const int index) const {
  assignVariable<T>(var_x, var_y, var_z, 1.0, keyword_query, sub_key, index);
}

} // namespace namelist
} // namespace stormm
