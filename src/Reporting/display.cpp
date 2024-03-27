#include <sys/ioctl.h>
#include <unistd.h>
#include "copyright.h"
#include "display.h"

namespace stormm {
namespace display {

//-------------------------------------------------------------------------------------------------
std::string horizontalRule(const std::string &left_corner, const std::string &right_corner,
                           const int width, const char middle) {
  const int n_left = left_corner.size();
  const int n_right = right_corner.size();
  const int actual_width = std::max(width, n_left + n_right);
  const int n_middle = actual_width - left_corner.size() - right_corner.size();
  std::string result(actual_width + 1, middle);
  for (int i = 0; i < n_left; i++) {
    result[i] = left_corner[i];
  }
  for (int i = 0; i < n_right; i++) {
    result[n_left + n_middle + i] = right_corner[i];
  }
  result[actual_width] = '\n';
  return result;
}

//-------------------------------------------------------------------------------------------------
void terminalHorizontalRule(const std::string &left_corner, const std::string &right_corner,
                            const int width, const char middle, std::ostream *foutp) {

  // Obtain the console size
  struct winsize console_dims;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &console_dims);
  int actual_width = (width == 0) ? console_dims.ws_col - 1 : width;
  const std::string hrule = horizontalRule(left_corner, right_corner, actual_width, middle);
  foutp->write(hrule.c_str(), hrule.size());
}

} // namespace display
} // namespace stormm
