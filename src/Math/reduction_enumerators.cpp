#include "copyright.h"
#include "reduction_enumerators.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ReductionStage input) {
  switch (input) {
  case ReductionStage::GATHER:
    return std::string("GATHER");
  case ReductionStage::SCATTER:
    return std::string("SCATTER");
  case ReductionStage::RESCALE:
    return std::string("RESCALE");
  case ReductionStage::ALL_REDUCE:
    return std::string("ALL_REDUCE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ReductionGoal input) {
  switch (input) {
  case ReductionGoal::NORMALIZE:
    return std::string("NORMALIZE");
  case ReductionGoal::CENTER_ON_ZERO:
    return std::string("CENTER_ON_ZERO");
  case ReductionGoal::CONJUGATE_GRADIENT:
    return std::string("CONJUGATE_GRADIENT");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const RdwuPerSystem input) {
  switch (input) {
  case RdwuPerSystem::ONE:
    return std::string("ONE");
  case RdwuPerSystem::MULTIPLE:
    return std::string("MULTIPLE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const RdwuAbstractMap input) {
  switch (input) {
  case RdwuAbstractMap::ATOM_START:
    return std::string("ATOM_START");
  case RdwuAbstractMap::ATOM_END:
    return std::string("ATOM_END");
  case RdwuAbstractMap::RESULT_INDEX:
    return std::string("RESULT_INDEX");    
  case RdwuAbstractMap::DEPN_START:
    return std::string("DEPN_START");    
  case RdwuAbstractMap::DEPN_END:
    return std::string("DEPN_END");    
  case RdwuAbstractMap::SYSTEM_ID:
    return std::string("SYSTEM_ID");    
  }
  __builtin_unreachable();
}

} // namespace stmath
} // namespace stormm
