#include "copyright.h"
#include "Reporting/error_format.h"
#include "textguard.h"

namespace stormm {
namespace parse {

//-------------------------------------------------------------------------------------------------
TextGuard::TextGuard(const std::string &left_in, const std::string &right_in,
                     const LineSpan spans_in) :
  left{left_in},
  right{right_in},
  left_length{static_cast<int>(left.size())},
  right_length{static_cast<int>(right.size())},
  line_span{spans_in},
  requires_termination{(right_in.size() > 0)}
{
  // Check to make sure that the initialization was done correctly
  if (line_span == LineSpan::MULTIPLE && right_length == 0) {
    rtErr("A text guard cannot apply to multiple lines without a right-hand delimiter.",
          "TextGuard", "TextGuard");
  }
  if (left_length == 0) {
    rtErr("A text guard must have a valid left-hand delimiter.");
  }
}

//-------------------------------------------------------------------------------------------------
const std::string& TextGuard::getLeft() const {
  return left;
}

//-------------------------------------------------------------------------------------------------
const std::string& TextGuard::getRight() const {
  return right;
}

//-------------------------------------------------------------------------------------------------
int TextGuard::leftSize() const {
  return left_length;
}

//-------------------------------------------------------------------------------------------------
int TextGuard::rightSize() const {
  return right_length;
}

//-------------------------------------------------------------------------------------------------
LineSpan TextGuard::getSpan() const {
  return line_span;
}

//-------------------------------------------------------------------------------------------------
bool TextGuard::getTerminationRequirement() const {
  return requires_termination;
}

} // namespace parse
} // namespace stormm
