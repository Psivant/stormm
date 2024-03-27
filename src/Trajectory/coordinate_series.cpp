#include "copyright.h"
#include "coordinate_series.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
ValidCoordinateTypes::ValidCoordinateTypes(const size_t double_id_in, const size_t float_id_in,
                                           const size_t short_id_in, const size_t int_id_in,
                                           const size_t llint_id_in) :
    double_id{double_id_in}, float_id{float_id_in}, short_id{short_id_in}, int_id{int_id_in},
    llint_id{llint_id_in}
{}

} // namespace trajectory
} // namespace stormm
