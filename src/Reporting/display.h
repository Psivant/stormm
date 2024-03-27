// -*-c++-*-
#ifndef STORMM_DISPLAY_H
#define STORMM_DISPLAY_H

#include <fstream>
#include <iostream>
#include <string>
#include "copyright.h"

namespace stormm {
namespace display {

/// \brief Produce a horizontal bar of repeating charaters using end points of the developer's
///        choice.
///
/// \param left_corner   The left border of the rule bar
/// \param right_corner  The right border of the fule bar
/// \param width         Total width of the rule bar, including the left and right end points
std::string horizontalRule(const std::string &left_corner = std::string("+"),
                           const std::string &right_corner = std::string("+"), int width = 79,
                           const char middle = '-');
  
/// \brief Print a horizontal bar of --- to a file (or the terminal screen), using end points of
///        the developer's choice (default '+')
///
/// \param left_corner   The left border of the rule bar
/// \param right_corner  The right border of the fule bar
/// \param width         Total width of the rule bar, including the left and right end points
/// \param foutp         The file stream (defaults to the terminal)
void terminalHorizontalRule(const std::string &left_corner = std::string("+"),
                            const std::string &right_corner = std::string("+"), int width = 0,
                            const char middle = '-', std::ostream *foutp = &std::cout);

} // namespace display
} // namespace stormm

#endif
