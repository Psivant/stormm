#include "copyright.h"
#include "coordinate_swap_plan.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
CoordinateSwapPlan::CoordinateSwapPlan(const CoordinateLineFormat layout_in,
                                       const int atom_count_in, const int total_width_in,
                                       const int fraction_in, const int x_start_in,
                                       const int y_start_in, const int z_start_in,
                                       const std::vector<int> &coordinate_lines_in,
                                       const std::vector<int> &x_widths_in,
                                       const std::vector<int> &y_widths_in,
                                       const std::vector<int> &z_widths_in,
                                       const std::vector<int> &x_starts_in,
                                       const std::vector<int> &y_starts_in,
                                       const std::vector<int> &z_starts_in) :
    layout{layout_in}, atom_count{atom_count_in}, total_width{total_width_in},
    fraction{fraction_in}, x_start{x_start_in}, y_start{y_start_in}, z_start{z_start_in},
    coordinate_lines{coordinate_lines_in}, x_widths{x_widths_in}, y_widths{y_widths_in},
    z_widths{z_widths_in}, x_starts{x_starts_in}, y_starts{y_starts_in}, z_starts{z_starts_in}
{}

} // namespace trajectory
} // namespace stormm
