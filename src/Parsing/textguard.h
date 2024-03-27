// -*-c++-*-
#ifndef STORMM_TEXTGUARD_H
#define STORMM_TEXTGUARD_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace parse {

/// \brief An enumerator for whether a text feature must take place on a single line or can span
///        multiple lines
enum class LineSpan {
  SINGLE, MULTIPLE
};

/// \brief Structure for specifying the features of a comment in some text file
class TextGuard {
public:

  /// \brief Constructor for text comment records
  ///
  /// \param left_in   Left-hand delimiter
  /// \param right_in  Right-hand delimiter
  /// \param spans_in  Line span rules (single or multiple lines)
  TextGuard(const std::string &left_in, const std::string &right_in = std::string(""),
            LineSpan spans_in = LineSpan::SINGLE);

  /// \brief Default destructor
  ~TextGuard() = default;

  /// \brief Get the left-hand guard
  const std::string& getLeft() const;

  /// \brief Get the right-hand guard
  const std::string& getRight() const;

  /// \brief Get the size of the left-hand guard, in terms of the number of characters
  int leftSize() const;

  /// \brief Get the size of the right-hand guard, in terms of the number of characters
  int rightSize() const;

  /// \brief Get the line span of this TextGuard object: is it limited to a single line, or can
  ///        it span multiple lines?
  LineSpan getSpan() const;

  /// \brief Get the (boolean) termination requirement: yes or no, does guarded text need to end
  ///        with an explicit right-hand guard? (If the TextGuard has no right-hand delimiter then
  ///        the answer is no, and the TextGuard must also be single-line in its scope.)
  bool getTerminationRequirement() const;

private:
  std::string left;          ///< The left hand delimiter (required)
  std::string right;         ///< The right hand delimiter (optional)
  int left_length;           ///< Character length of the left delimiter
  int right_length;          ///< Character length of the right delimiter
  LineSpan line_span;        ///< May be MULTIPLE if "right" is not blank, although some
                             ///<   left / right delimiter pairs will have SINGLE spans
  bool requires_termination; ///< True if the string "right" is not blank
};

} // namespace parse
} // namespace stormm

#endif
