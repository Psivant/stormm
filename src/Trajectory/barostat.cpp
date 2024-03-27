#include "copyright.h"
#include "barostat.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
Barostat::Barostat() :
    kind{BarostatKind::NONE}
{}

//-------------------------------------------------------------------------------------------------
Barostat::Barostat(const BarostatKind kind_in) :
    kind{kind_in}
{}

//-------------------------------------------------------------------------------------------------
BarostatKind Barostat::getKind() const {
  return kind;
}
  
//-------------------------------------------------------------------------------------------------
std::string getBarostatName(const BarostatKind kind) {
  switch (kind) {
  case BarostatKind::NONE:
    return std::string("NONE");
  case BarostatKind::MONTE_CARLO:
    return std::string("MONTE_CARLO");
  case BarostatKind::BERENDSEN:
    return std::string("BERENDSEN");
  }
  __builtin_unreachable();
}

} // namespace trajectory
} // namespace stormm
