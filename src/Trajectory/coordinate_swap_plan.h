// -*-c++-*-
#ifndef STORMM_COORDINATE_SWAP_PLAN_H
#define STORMM_COORDINATE_SWAP_PLAN_H

#include <vector>
#include "copyright.h"
#include "trajectory_enumerators.h"

namespace stormm {
namespace trajectory {

/// \brief Collect the information needed to swap the coordinates (and, if applicable, the box
///        information) in a highly annotated block of ascii text.
class CoordinateSwapPlan {
public:

  /// \brief The constructors can accept any and all inputs.
  CoordinateSwapPlan(CoordinateLineFormat layout = CoordinateLineFormat::FIXED_COLUMN,
                     int atom_count_in = 0, int total_width_in = 10, int fraction_in = 4,
                     int x_start_in = 0, int y_start_in = 10, int z_start_in = 20,
                     const std::vector<int> &coordinate_lines_in = {},
                     const std::vector<int> &x_widths_in = {},
                     const std::vector<int> &y_widths_in = {},
                     const std::vector<int> &z_widths_in = {},
                     const std::vector<int> &x_starts_in = {},
                     const std::vector<int> &y_starts_in = {},
                     const std::vector<int> &z_starts_in = {});
  
private:
  CoordinateLineFormat layout;        ///< Layout of the text: fixed column or free format
  int atom_count;                     ///< Number of atoms, and the length of arrays in this object
  int total_width;                    ///< Total width of each number in a fixed-column format
  int fraction;                       ///< Number of digits after the decimal
  int x_start;                        ///< General starting location for X coordinates, in a
                                      ///<   fixed-column format
  int y_start;                        ///< General starting location for Y coordinates
  int z_start;                        ///< General starting location for Z coordinates
  std::vector<int> coordinate_lines;  ///< A list of integers holding coordinate lines in the
                                      ///<   associated text file (indexing starts at zero)
  std::vector<int> x_widths;          ///< Widths of X coordinate numbers (if the lines do not
                                      ///<   keep to a fixed column format)
  std::vector<int> y_widths;          ///< Widths of Y coordinate numbers
  std::vector<int> z_widths;          ///< Widths of Z coordinate numbers
  std::vector<int> x_starts;          ///< Starting character indices for X coordinate numbers (if
                                      ///<   the lines do not keep to a fixed column format)
  std::vector<int> y_starts;          ///< Starting character indices for Y coordinate numbers
  std::vector<int> z_starts;          ///< Starting character indices for Z coordinate numbers
};
  
} // namespace trajectory
} // namespace stormm

#endif
