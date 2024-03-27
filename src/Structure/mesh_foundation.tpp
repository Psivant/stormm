// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename T>
const CoordinateSeries<T>* MeshFoundation::getEnsemblePointer() const {
  return reinterpret_cast<CoordinateSeries<T>*>(cf_ensemble);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MeshFoundation::setEnsemble(const CoordinateSeries<T> *cs_in) {
  const CoordinateSeries<void> *cs_voided = reinterpret_cast<const CoordinateSeries<void>*>(cs_in);
  cf_ensemble = const_cast<CoordinateSeries<void>*>(cs_voided);
  ensemble_data_type = std::type_index(typeid(T)).hash_code();
  testReadiness();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MeshFoundation::setEnsemble(const CoordinateSeries<T> &cs_in) {
  setEnsemble(cs_in.getSelfPointer());
  testReadiness();
}

} // namespace structure
} // namespace stormm
