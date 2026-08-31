#include "copyright.h"
#include "debugging_enumerators.h"

namespace stormm {
namespace debug {

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const AnomalyContext input) {
  switch (input) {
  case AnomalyContext::PHASESPACE_SYNTHESIS:
    return std::string("PHASESPACE_SYNTHESIS");
  case AnomalyContext::CELLGRID:
    return std::string("CELLGRID");
  case AnomalyContext::QQ_CELLGRID:
    return std::string("QQ_CELLGRID");
  case AnomalyContext::LJ_CELLGRID:
    return std::string("LJ_CELLGRID");
  }
  __builtin_unreachable();
}

} // namespace debug
} // namespace stormm
