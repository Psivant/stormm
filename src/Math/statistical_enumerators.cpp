#include "statistical_enumerators.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const VarianceMethod input) {
  switch (input) {
  case VarianceMethod::VARIANCE:
    return std::string("VARIANCE");
  case VarianceMethod::STANDARD_DEVIATION:
    return std::string("STANDARD_DEVIATION");
  case VarianceMethod::ROOT_MEAN_SQUARED_DEVIATION:
    return std::string("ROOT_MEAN_SQUARED_DEVIATION");
  case VarianceMethod::COEFFICIENT_OF_VARIATION:
    return std::string("COEFFICIENT_OF_VARIATION");
  case VarianceMethod::NORMALIZED_RMSD:
    return std::string("NORMALIZED_RMSD");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const DataOrder input) {
  switch (input) {
  case DataOrder::ASCENDING:
    return std::string("ASCENDING");
  case DataOrder::DESCENDING:
    return std::string("DESCENDING");
  case DataOrder::NONE:
    return std::string("NONE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ResamplingMethod input) {
  switch (input) {
  case ResamplingMethod::JACKKNIFE:
    return std::string("JACKKNIFE");
  case ResamplingMethod::BOOTSTRAP:
    return std::string("BOOTSTRAP");
  }
  __builtin_unreachable();
}

} // namespace stmath
} // namespace stormm
