// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace analysis {

//-------------------------------------------------------------------------------------------------
template <typename T>
ComparisonGuide::ComparisonGuide(const CoordinateSeries<T> *cs_in, const AtomGraph *ag_in,
                                 const GpuDetails &gpu) :
    ComparisonGuide()
{
  cs_ptr = reinterpret_cast<CoordinateSeries<int>*>(const_cast<CoordinateSeries<T>*>(cs_in));
  basis = StructureSource::SERIES;
  csptr_data_type = std::type_index(typeid(T)).hash_code();
  ag_pointers.resize(1);
  ag_pointers[0] = const_cast<AtomGraph*>(ag_in);
  allocateSystemIndexing();
  setSystemIndexing();
  setWorkUnits(gpu);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
ComparisonGuide::ComparisonGuide(const CoordinateSeries<T> &cs_in, const AtomGraph &ag_in,
                                 const GpuDetails &gpu) :
    ComparisonGuide(cs_in.getSelfPointer(), ag_in.getSelfPointer(), gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
const CoordinateSeries<T>* ComparisonGuide::getSeriesPointer() const {
  if (cs_ptr == nullptr) {
    rtErr("No CoordinateSeries is referenced.", "ComparsionGuide", "getSeriesPointer");
  }
  return reinterpret_cast<CoordinateSeries<T>*>(cs_ptr);
}

} // namespace analysis
} // namespace stormm
