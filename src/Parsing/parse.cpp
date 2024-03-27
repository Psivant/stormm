#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>
#include "copyright.h"
#include "Reporting/error_format.h"
#include "ascii_numbers.h"
#include "parse.h"

namespace stormm {
namespace parse {

//-------------------------------------------------------------------------------------------------
std::string operator+(const std::string &lhs, const char4 rhs) {
  std::string result = lhs;
  if (rhs.x == '\0') return result;
  result += rhs.x;
  if (rhs.y == '\0') return result;
  result += rhs.y;
  if (rhs.z == '\0') return result;
  result += rhs.z;
  if (rhs.w == '\0') return result;
  result += rhs.w;
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string char2ToString(const char2 value) {
  std::string result;
  if (value.x == '\0') return result;
  result += value.x;
  if (value.y == '\0') return result;
  result += value.y;
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string char3ToString(const char3 value) {
  std::string result;
  if (value.x == '\0') return result;
  result += value.x;
  if (value.y == '\0') return result;
  result += value.y;
  if (value.z == '\0') return result;
  result += value.z;
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string char4ToString(const char4 value) {
  std::string result;
  if (value.x == '\0') return result;
  result += value.x;
  if (value.y == '\0') return result;
  result += value.y;
  if (value.z == '\0') return result;
  result += value.z;
  if (value.w == '\0') return result;
  result += value.w;
  return result;
}

//-------------------------------------------------------------------------------------------------
char4 stringToChar4(const std::string &value) {
  char4 result;
  const int vlen = value.size();
  result.x = (vlen > 0) ? value[0] : ' ';
  result.y = (vlen > 1) ? value[1] : ' ';
  result.z = (vlen > 2) ? value[2] : ' ';
  result.w = (vlen > 3) ? value[3] : ' ';
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string alphabetNumber(ullint input) {
  int nchar = 1;
  ullint factor = 26;
  ullint range = factor;
  while (range < input) {
    factor *= 26;
    range += factor;
    nchar++;
  }
  std::string result(nchar, ' ');
  ullint remainder = input;
  factor /= 26;
  for (int i = nchar - 1; i >= 0; i--) {
    const ullint digit = remainder / factor;
    remainder -= digit * factor;
    result[nchar - 1 - i] = 'a' + digit - (i > 0);
    factor /= 26;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string rgbHexCode(const uchar4 color_in) {
  std::string result(6, '0');
  const std::vector<int> hexdigits = {
    static_cast<int>((color_in.x >> 4) & 0xf), static_cast<int>(color_in.x & 0xf),
    static_cast<int>((color_in.y >> 4) & 0xf), static_cast<int>(color_in.y & 0xf),
    static_cast<int>((color_in.z >> 4) & 0xf), static_cast<int>(color_in.z & 0xf)
  };
  const int ascii_zero = '0';
  const int ascii_nine = '9';
  const int ascii_a = 'a';
  for (size_t i = 0; i < 6; i++) {
    result[i] = (hexdigits[i] >= 0 && hexdigits[i] <= 9) ? hexdigits[i] + ascii_zero :
                                                           hexdigits[i] - 10 + ascii_a;
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
bool verifyNumberFormat(const char* a, const NumberFormat cform, const int read_begin,
                        const int len) {

  // Pre-allocate a buffer and fill it with string termination characters
  const int width = (len > 0) ? len : strlen(a) - read_begin;
  std::vector<char> buffer(width, ' ');

  // The CHAR4 type can accept any characters
  if (cform == NumberFormat::CHAR4) {
    return true;
  }
  
  // Fill the buffer with data from the input character array
  bool problem = false;
  bool number_begins = false;
  bool number_ends = false;
  for (int i = 0; i < width; i++) {
    buffer[i] = a[read_begin + i];
    number_begins = (number_begins || buffer[i] != ' ');
    number_ends = (number_ends || (number_begins && buffer[i] == ' '));
    problem = (problem || (buffer[i] != ' ' && number_ends));
  }

  // Check each character
  bool e_found = false;
  bool dot_found = false;
  int signs_found = 0;
  int minus_found = 0;
  int digit_found = 0;
  for (int i = 0; i < width; i++) {
    if (buffer[i] == 'E' || buffer[i] == 'e') {
      problem = (problem || e_found);
      e_found = true;
    }
    else if (buffer[i] == '.') {
      problem = (problem || dot_found || e_found);
      dot_found = true;
    }
    else if (buffer[i] == '+' || buffer[i] == '-') {
      problem = (problem || (i > 0 && buffer[i - 1] != ' ' && e_found == false));
      minus_found += (buffer[i] == '-');
      signs_found++;
    }
    else {
      switch(cform) {
      case NumberFormat::SCIENTIFIC:
      case NumberFormat::STANDARD_REAL:
      case NumberFormat::INTEGER:
      case NumberFormat::LONG_LONG_INTEGER:
      case NumberFormat::UNSIGNED_INTEGER:
      case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
        problem = (problem || ((buffer[i] < '0' || buffer[i] > '9') && buffer[i] != ' '));
        break;
      case NumberFormat::CHAR4:
        break;
      }
    }
  }
  problem = (problem || number_begins == false);

  // Note any malformed numbers
  switch(cform) {
  case NumberFormat::SCIENTIFIC:
    problem = (problem || (e_found == false || dot_found == false || signs_found > 2));
    break;
  case NumberFormat::STANDARD_REAL:
    problem = (problem || (signs_found > 1 || e_found));
    break;
  case NumberFormat::INTEGER:
  case NumberFormat::LONG_LONG_INTEGER:
    problem = (problem || (signs_found > 1 || dot_found || e_found));
    break;
  case NumberFormat::UNSIGNED_INTEGER:
  case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
    problem = (problem || (signs_found > 1 || minus_found > 0 || dot_found || e_found));
    break;
  case NumberFormat::CHAR4:
    problem = (problem || width != 4);
    break;
  }

  return (problem == false);
}

//-------------------------------------------------------------------------------------------------
bool verifyContents(const char* line, const int start_pos, const int length,
                    const NumberFormat fmt) {
  return verifyNumberFormat(line, fmt, start_pos, length);
}

//-------------------------------------------------------------------------------------------------
bool verifyContents(const std::string &line, const NumberFormat fmt) {
  return verifyNumberFormat(line.c_str(), fmt, 0, line.size());
}

//-------------------------------------------------------------------------------------------------
bool verifyContents(const std::string &line, const int start_pos, const int length,
                    const NumberFormat fmt) {
  const int line_length = line.size();
  if (start_pos >= line_length || start_pos + length >= line_length) {
    return false;
  }
  else {
    return verifyNumberFormat(line.c_str(), fmt, start_pos, length);
  }
}

//-------------------------------------------------------------------------------------------------
bool verifyContents(const TextFile &tf, const int line, const int start_pos, const int length,
                    const NumberFormat fmt) {
  const int line_length = tf.getLineLength(line);
  if (start_pos >= line_length || start_pos + length >= line_length) {
    return false;
  }
  else {
    return verifyNumberFormat(tf.getLinePointer(line), fmt, start_pos, length);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool verifyContents(const TextFileReader &tfr, const int line, const int start_pos,
                    const int length, const NumberFormat fmt) {
  if (line >= tfr.line_count || start_pos >= tfr.line_lengths[line] ||
      start_pos + length >= tfr.line_lengths[line]) {
    return false;
  }
  else {
    return verifyNumberFormat(&tfr.text[tfr.line_limits[line]], fmt, start_pos, length);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
char uppercase(const char tc) {
  return tc - (tc >= 97 && tc <= 122) * 32;
}

//-------------------------------------------------------------------------------------------------
std::string uppercase(const std::string &ts) {
  std::string result = ts;
  const size_t n_char = ts.size();
  for (size_t i = 0; i < n_char; i++) {
    const char tc = result[i];
    result[i] = tc - (tc >= 97 && tc <= 122) * 32;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void uppercase(char* tcs, const size_t n_char) {
  const size_t scan_char = (n_char == 0) ? strlen(tcs) : n_char;
  for (size_t i = 0; i < scan_char; i++) {
    const char tc = tcs[i];
    tcs[i] = tc - (tc >= 97 && tc <= 122) * 32;
  }
}

//-------------------------------------------------------------------------------------------------
std::string uppercase(const char* tcs) {
  std::string result(tcs);
  const size_t n_char = result.size();
  for (size_t i = 0; i < n_char; i++) {
    const char tc = tcs[i];
    result[i] = tc - (tc >= 97 && tc <= 122) * 32;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
char lowercase(const char tc) {
  return tc + (tc >= 65 && tc <= 90) * 32;
}

//-------------------------------------------------------------------------------------------------
std::string lowercase(const std::string &ts) {
  std::string result = ts;
  const size_t n_char = ts.size();
  for (size_t i = 0; i < n_char; i++) {
    const char tc = result[i];
    result[i] = tc + (tc >= 65 && tc <= 90) * 32;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void lowercase(char* tcs, const size_t n_char) {
  const size_t scan_char = (n_char == 0) ? strlen(tcs) : n_char;
  for (size_t i = 0; i < scan_char; i++) {
    const char tc = tcs[i];
    tcs[i] = tc + (tc >= 65 && tc <= 90) * 32;
  }
}

//-------------------------------------------------------------------------------------------------
std::string lowercase(const char* tcs) {
  std::string result(tcs);
  const size_t n_char = result.size();
  for (size_t i = 0; i < n_char; i++) {
    const char tc = tcs[i];
    result[i] = tc + (tc >= 65 && tc <= 90) * 32;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
bool strcmpCased(const char* sa, const char* sb, const CaseSensitivity csen) {
  int i = 0;
  switch (csen) {
  case CaseSensitivity::YES:
    while (sa[i] == sb[i]) {
      if (sa[i] == '\0') {
        return true;
      }
      i++;
    }
    return false;
  case CaseSensitivity::NO:
    while (uppercase(sa[i]) == uppercase(sb[i])) {
      if (sa[i] == '\0') {
        return true;
      }
      i++;
    }
    return false;
  case CaseSensitivity::AUTOMATIC:
    rtErr("No AUTOMATIC behavior is defined for case-based string comparison.  AUTOMATIC "
          "settings for case sensitivity are defined at higher levels for specific situations, "
          "not the low-level implementation.", "strcmpCased");
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool strcmpCased(const std::string &sa, const char* sb, const CaseSensitivity csen) {
  return strcmpCased(sa.c_str(), sb, csen);
}

//-------------------------------------------------------------------------------------------------
bool strcmpCased(const char* sa, const std::string &sb, const CaseSensitivity csen) {
  return strcmpCased(sa, sb.c_str(), csen);
}

//-------------------------------------------------------------------------------------------------
bool strcmpCased(const std::string &sa, const std::string &sb, const CaseSensitivity csen) {
  return strcmpCased(sa.c_str(), sb.c_str(), csen);
}

//-------------------------------------------------------------------------------------------------
bool strncmpCased(const char* sa, const char* sb, const int length, const CaseSensitivity csen) {
  bool match = true;
  switch (csen) {
  case CaseSensitivity::YES:
    for (int i = 0; i < length; i++) {
      match = (match && (sa[i] == sb[i]));
      if (sa[i] == '\0' && match) {
        return true;
      }
    }
    break;
  case CaseSensitivity::NO:
    for (int i = 0; i < length; i++) {
      match = (match && (uppercase(sa[i]) == uppercase(sb[i])));
      if (sa[i] == '\0' && match) {
        return true;
      }
    }
    break;
  case CaseSensitivity::AUTOMATIC:
    rtErr("No AUTOMATIC behavior is defined for case-based string comparison.  AUTOMATIC "
          "settings for case sensitivity are defined at higher levels for specific situations, "
          "not the low-level implementation.", "strncmpCased");
    break;
  }
  return match;
}

//-------------------------------------------------------------------------------------------------
bool strncmpCased(const std::string &sa, const char* sb, const CaseSensitivity csen,
                  const int length) {
  const int actual_length = (length < 0) ? sa.size() : length;
  if (length > static_cast<int>(sa.size())) {
    return false;
  }
  return strncmpCased(sa.c_str(), sb, actual_length, csen);
}

//-------------------------------------------------------------------------------------------------
bool strncmpCased(const char* sa, const std::string &sb, const CaseSensitivity csen,
                  const int length) {
  const int actual_length = (length < 0) ? sb.size() : length;
  if (length > static_cast<int>(sb.size())) {
    return false;
  }
  return strncmpCased(sa, sb.c_str(), actual_length, csen);
}

//-------------------------------------------------------------------------------------------------
bool strncmpCased(const std::string &sa, const std::string &sb, const CaseSensitivity csen,
                  const int length) {
  if (length < 0 && sa.size() != sb.size()) {
    rtErr("The extent of a comparison cannot be inferred from two string arguments of different "
          "lengths (" + std::to_string(sa.size()) + " and " + std::to_string(sb.size()) + ", " +
          sa + " and " + sb + ").", "strncmpCased");
  }
  if (length > static_cast<int>(sa.size()) || length > static_cast<int>(sb.size())) {
    rtErr("String of length " + std::to_string(sa.size()) + " and " + std::to_string(sb.size()) +
          " cannot be subjected to a comparison of " + std::to_string(length) + " characters.",
          "strncmpCased");
  }
  const int actual_length = (length < 0) ? sa.size() : length;
  return strncmpCased(sa.c_str(), sb.c_str(), actual_length, csen);
}

//-------------------------------------------------------------------------------------------------
bool strcmpWildCard(const std::string &target, const std::string &query,
                    const std::vector<WildCardKind> &wildcards) {

  // Check that the wildcards array is of the appropriate length
  if (wildcards.size() != query.size()) {
    rtErr("A wildcard string comparison requires that the second string, with " +
          std::to_string(query.size()) + " characters, be supported by the same number of "
          "wildcards (currently " + std::to_string(wildcards.size()) + ").", "strncmpWildCard");
  }
  
  // Ignore white space at the beginning of the target
  int target_con = 0;
  const int target_len = target.size();
  while (target_con < target_len && target[target_con] == ' ') {
    target_con++;
  }
  int query_con = 0;
  const int query_len = query.size();
  bool star_power_buff = false;
  const int n_wild = wildcards.size();
  for (int j = 0; j < query_len; j++) {

    // White space in the query does nothing to advance matching in the target
    if (query[j] == ' ') {
      query_con++;
      continue;
    }
    switch (wildcards[j]) {
    case WildCardKind::NONE:
      if (star_power_buff) {
        while (target_con < target_len && query[j] != target[target_con]) {
          target_con++;
        }
        star_power_buff = false;
      }
      if (target_con < target_len && query[j] == target[target_con]) {
        query_con++;
        target_con++;
      }
      break;
    case WildCardKind::FREE_CHARACTER:
      if (target_con < target_len) {
        query_con++;
        target_con++;
      }
      break;
    case WildCardKind::FREE_STRETCH:
      query_con++;
      star_power_buff = true;
      break;
    }
  }

  // If the query ends with a free stretch wildcard active, keep it going until the end
  target_con = (star_power_buff) ? target_len : target_con;

  // Ignore white space at the end of the target
  while (target_con < target_len && target[target_con] == ' ') {
    target_con++;
  }
  
  return (target_con == target_len && query_con == query_len);
}

//-------------------------------------------------------------------------------------------------
std::string addLeadingWhiteSpace(const std::string &input, const size_t target_length) {
  if (input.size() >= target_length) {
    return input;
  }
  std::string result;
  result.reserve(target_length + target_length - input.size());
  const size_t nadd = target_length - input.size();
  for (size_t i = 0; i < nadd; i++) {
    result.append(" ");
  }
  result.append(input);
  return result;
}

//-------------------------------------------------------------------------------------------------
void addLeadingWhiteSpace(std::string *input, const size_t target_length) {
  if (input->size() >= target_length) {
    return;
  }
  const size_t init_length = input->size();
  const size_t nadd = target_length - init_length;
  input->resize(target_length);
  char* iptr = input->data();
  for (size_t i = init_length - 1; i < init_length; i--) {
    iptr[i + nadd] = iptr[i];
  }
  for (size_t i = 0; i < nadd; i++) {
    iptr[i] = ' ';
  }
}

//-------------------------------------------------------------------------------------------------
std::string addTailingWhiteSpace(const std::string &input, const size_t target_length) {
  if (input.size() >= target_length) {
    return input;
  }
  std::string result;
  result.reserve(target_length);
  result.append(input);
  const size_t nadd = target_length - input.size();
  for (size_t i = 0; i < nadd; i++) {
    result.append(" ");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void addTailingWhiteSpace(std::string *input, const size_t target_length) {
  if (input->size() >= target_length) {
    return;
  }
  const size_t init_length = input->size();
  input->resize(target_length);
  char* iptr = input->data();
  for (size_t i = init_length; i < target_length; i++) {
    iptr[i] = ' ';    
  }
}

//-------------------------------------------------------------------------------------------------
std::string removeLeadingWhiteSpace(const std::string &input) {
  size_t leading_space = 0;
  const size_t slen = input.size();
  while (leading_space < slen && input[leading_space] == ' ') {
    leading_space++;
  }
  if (leading_space > 0) {
    return input.substr(leading_space);
  }
  else {
    return input;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void removeLeadingWhiteSpace(std::string *input) {
  size_t leading_space = 0;
  const size_t slen = input->size();
  char* iptr = input->data();
  while (leading_space < slen && iptr[leading_space] == ' ') {
    leading_space++;
  }
  if (leading_space > 0) {
    for (size_t i = leading_space; i < slen; i++) {
      iptr[i - leading_space] = iptr[i];
    }
    input->resize(slen - leading_space);
  }
}

//-------------------------------------------------------------------------------------------------
std::string removeTailingWhiteSpace(const std::string &input) {
  const size_t slen = input.size();
  size_t tailing_space = slen - 1;
  while (tailing_space < slen && input[tailing_space] == ' ') {
    tailing_space--;
  }
  if (tailing_space < slen) {
    return input.substr(0, tailing_space + 1);
  }
  else {
    return std::string("");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void removeTailingWhiteSpace(std::string *input) {
  const size_t slen = input->size();
  size_t tailing_space = slen - 1;
  const char* iptr = input->data();
  while (tailing_space < slen && iptr[tailing_space] == ' ') {
    tailing_space--;
  }
  if (tailing_space < slen) {
    input->resize(tailing_space + 1);
  }
  else {
    input->resize(0);
  }
}

//-------------------------------------------------------------------------------------------------
size_t justifyStrings(std::vector<std::string> *ds, const int lower_limit,
                      const int upper_limit, const JustifyText fside, const size_t max_length,
                      const int long_count, const int near_limit) {

  // Some parameters may have very long defaults, or multiple values that result in a very long
  // string.  Catch those cases and go with a consensus length that suits most parameters to
  // avoid polluting the output with excessive white space.
  size_t dflt_length = 0;
  std::string* ds_ptr = ds->data();
  if (lower_limit >= upper_limit || lower_limit > ds->size() || upper_limit > ds->size()) {
    rtErr("The range " + std::to_string(lower_limit) + " - " + std::to_string(upper_limit) +
          " is invalid for a list of " + std::to_string(ds->size()) + " strings.",
          "justifyStrings");
  }
  for (int i = lower_limit; i < upper_limit; i++) {
    dflt_length = std::max(dflt_length, ds_ptr[i].size());
  }
  if (dflt_length >= max_length) {

    // Check to see whether there are a sufficient number of strings with lengths near that of
    // the longest string.  If not, clip the longest strings, up to long_count, until there are
    // a sufficient number of strings with nearly the same length.
    const int n_params = upper_limit - lower_limit;
    std::vector<size_t> all_lengths(n_params);
    for (size_t i = lower_limit; i < upper_limit; i++) {
      all_lengths[i - lower_limit] = ds_ptr[i].size();
    }
    std::sort(all_lengths.begin(), all_lengths.end(), [](int x, int y) { return x < y; });
    int pivot = n_params - 1;
    bool critical_mass;
    do {
      const int llim = std::min(std::max(pivot - long_count + 1, 0), pivot);
      if (pivot < n_params - long_count || all_lengths[llim] + near_limit >= all_lengths[pivot]) {
        critical_mass = true;
      }
      else {
        pivot--;
        critical_mass = false;
      }
    } while (critical_mass == false);
    dflt_length = all_lengths[pivot];
  }
  for (int i = lower_limit; i < upper_limit; i++) {
    ds_ptr[i] = removeLeadingWhiteSpace(ds_ptr[i]);
    ds_ptr[i] = removeTailingWhiteSpace(ds_ptr[i]);
  }
  switch (fside) {
  case JustifyText::LEFT:
    for (int i = lower_limit; i < upper_limit; i++) {
      ds_ptr[i] = addTailingWhiteSpace(ds_ptr[i], dflt_length);
    }
    break;
  case JustifyText::RIGHT:
    for (int i = lower_limit; i < upper_limit; i++) {
      ds_ptr[i] = addLeadingWhiteSpace(ds_ptr[i], dflt_length);
    }
    break;
  case JustifyText::CENTER:
    for (int i = lower_limit; i < upper_limit; i++) {
      const int nspc = dflt_length - static_cast<int>(ds_ptr[i].size());
      if (nspc > 0) {
        const int lspc = nspc / 2;
        const int rspc = nspc - lspc;
        ds_ptr[i] = addTailingWhiteSpace(addLeadingWhiteSpace(ds_ptr[i], dflt_length - rspc),
                                         dflt_length);
      }
    }
    break;
  }
  return dflt_length;
}

//-------------------------------------------------------------------------------------------------
size_t justifyStrings(std::vector<std::string> *ds, const JustifyText fside,
                      const size_t max_length, const int long_count, const int near_limit) {
  return justifyStrings(ds, 0, ds->size(), fside, max_length, long_count, near_limit);
}

//-------------------------------------------------------------------------------------------------
size_t maximumLineLength(const std::string &input) {
  size_t llim = 0;
  size_t lmax = 0;
  const size_t slen = input.size();
  for (size_t i = 0; i < slen; i++) {
    if (input[i] == '\n') {
      lmax = std::max(lmax, i - llim);
      llim = i;
    }
  }
  lmax = std::max(lmax, slen - llim);
  return lmax;
}
//-------------------------------------------------------------------------------------------------
int realDecimalPlaces(const double value, const int limit) {
  
  // Determine the fractional component
  double abs_val = fabs(value);
  double frac = abs_val - floor(abs_val);
  double smallest_significant_amount = pow(0.1, limit);
  int n_place = 0;
  while (n_place < limit && 1.01 * frac > smallest_significant_amount &&
         1.01 * fabs(1.0 - frac) > smallest_significant_amount) {
    frac *= 10.0;
    frac -= floor(frac);
    n_place++;
  }
  return n_place;
}

//-------------------------------------------------------------------------------------------------
std::string realToString(const double value, const int format_a, const int format_b,
                         const NumberFormat method, const NumberPrintStyle style) {

  // Check the overall format
  switch (method) {
  case NumberFormat::SCIENTIFIC:
    if (format_a != free_number_format && format_b != free_number_format &&
        format_b > format_a - 7) {
      rtErr("A real number of format %" + std::to_string(format_a) + "." +
            std::to_string(format_b) + "e has too many decimal places for its overall length.",
            "realToString");
    }
    break;
  case NumberFormat::STANDARD_REAL:
    if (format_a != free_number_format && format_b != free_number_format &&
        format_b > format_a - 2) {
      rtErr("A real number of format %" + std::to_string(format_a) + "." +
            std::to_string(format_b) + "f has too many decimal places for its overall length.",
            "realToString");
    }
    break;
  case NumberFormat::INTEGER:
  case NumberFormat::LONG_LONG_INTEGER:
  case NumberFormat::UNSIGNED_INTEGER:
  case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
  case NumberFormat::CHAR4:
    rtErr("The printing method for a real number must be either SCIENTIFIC or STANDARD_REAL.",
          "realToString");
  }

  // Examine the given dimensions
  if (format_b < 0 && format_b != free_number_format) {
    rtErr("A nonsensical number of decimal places (" + std::to_string(format_b) +
          ") was specified.", "realToString");
  }
  if (abs(format_a) >= 62 && format_a != free_number_format) {
    rtErr("The requested number exceeds format limits (maximum 64 characters, %" +
          std::to_string(format_a) + "." + std::to_string(format_b) + "f requested).",
          "realToString");
  }

  // Print the formatted number and return it as a string
  char buffer[64];
  if (format_a != free_number_format && format_b != free_number_format) {
    switch (style) {
    case NumberPrintStyle::STANDARD:
      if (method == NumberFormat::SCIENTIFIC) {
        snprintf(buffer, 64, "%*.*e", format_a, format_b, value);
      }
      else {
        snprintf(buffer, 64, "%*.*f", format_a, format_b, value);
      }
      break;
    case NumberPrintStyle::LEADING_ZEROS:
      if (method == NumberFormat::SCIENTIFIC) {
        snprintf(buffer, 64, "%0*.*e", format_a, format_b, value);
      }
      else {
        snprintf(buffer, 64, "%0*.*f", format_a, format_b, value);
      }
      break;
    }
  }
  else if (format_b != free_number_format) {
    switch (style) {
    case NumberPrintStyle::STANDARD:
      if (method == NumberFormat::SCIENTIFIC) {
        snprintf(buffer, 64, "%.*e", format_b, value);
      }
      else {
        snprintf(buffer, 64, "%.*f", format_b, value);
      }
      break;
    case NumberPrintStyle::LEADING_ZEROS:
      if (method == NumberFormat::SCIENTIFIC) {
        snprintf(buffer, 64, "%0.*e", format_b, value);
      }
      else {
        snprintf(buffer, 64, "%0.*f", format_b, value);
      }
      break;
    }
  }
  else if (format_a != free_number_format) {
    switch (style) {
    case NumberPrintStyle::STANDARD:
      if (method == NumberFormat::SCIENTIFIC) {
        snprintf(buffer, 64, "%.*e", format_a, value);
      }
      else {
        snprintf(buffer, 64, "%.*f", format_a, value);

      }
      break;
    case NumberPrintStyle::LEADING_ZEROS:
      if (method == NumberFormat::SCIENTIFIC) {
        snprintf(buffer, 64, "%0.*e", format_a, value);
      }
      else {
        snprintf(buffer, 64, "%0.*f", format_a, value);
      }
      break;
    }
  }
  else {
    snprintf(buffer, 64, "%.*f", realDecimalPlaces(value, 8), value);
  }
  return std::string(buffer);
}

//-------------------------------------------------------------------------------------------------
std::string realToString(const double value, const int format_a, const NumberFormat method,
                         const NumberPrintStyle style) {
  return realToString(value, format_a, free_number_format, method, style);
}

//-------------------------------------------------------------------------------------------------
std::string minimalRealFormat(const double value, const double rel) {

  // Return immediately if the value is zero
  if (value == 0.0) {
    return std::string("0.0");
  }

  // Compute the critical increment
  const double increment = fabs(value * rel);
  
  // Find the best result possible in scientific notation
  int nsci_dec = ceil(fabs(log10(fabs(value)) - log10(fabs(value * rel)))) + 1;
  std::string best_sci_result;
  bool trim_success = false;
  do {
    const std::string sci_result = realToString(value, nsci_dec + 7, nsci_dec,
                                                NumberFormat::SCIENTIFIC);
    const double test = std::stod(sci_result);
    const double test_margin = fabs(value - test);
    if (test_margin <= increment) {
      if (nsci_dec > 1) {
        nsci_dec--;
        trim_success = true;
      }
      else {
        trim_success = false;
      }
      best_sci_result = sci_result;
    }
    else {
      trim_success = false;
    }
  } while (nsci_dec >= 1 && trim_success);
  best_sci_result = removeLeadingWhiteSpace(best_sci_result);
  const int nsci_char = best_sci_result.size();
  int nelim;
  for (int i = 0; i < nsci_char; i++) {
    if (best_sci_result[i] == 'e' && i < nsci_char - 2) {
      if (best_sci_result[i + 1] == '+' && best_sci_result[i + 2] == '0') {
        for (int j = i + 1; j < nsci_char - 2; j++) {
          best_sci_result[j] = best_sci_result[j + 2];
        }
        nelim = 2;
      }
      else if (best_sci_result[i + 2] == '0') {
        for (int j = i + 2; j < nsci_char - 1; j++) {
          best_sci_result[j] = best_sci_result[j + 1];
        }
        nelim = 1;
      }
    }
  }
  best_sci_result.resize(best_sci_result.size() - nelim);

  // Find the best result possible in standard notation
  int nstd_dec = 1;
  int nstd_maj = ceil(fabs(log10(fabs(value))));
  std::string best_std_result;
  bool search_on = true;
  while (search_on) {
    best_std_result = realToString(value, 2 + nstd_maj + nstd_dec, nstd_dec,
                                   NumberFormat::STANDARD_REAL);
    const double test = std::stod(best_std_result);
    const double test_margin = fabs(value - test);
    if (test_margin > increment) {
      nstd_dec++;
      search_on = true;
    }
    else {
      search_on = false;
    }
  }
  const size_t std_test_size = best_std_result.size();
  if (std_test_size > 2 && best_std_result[std_test_size - 1] == '0' &&
      best_std_result[std_test_size - 2] == '.') {
    best_std_result.resize(std_test_size - 2);
  }
  best_std_result = removeLeadingWhiteSpace(best_std_result);

  // Return the shorter representation, preferring standard notation in the case of a tie.
  if (best_std_result.size() <= best_sci_result.size()) {
    return best_std_result;
  }
  else {
    return best_sci_result;
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
std:: string intToString(const llint value, const int width, const NumberPrintStyle style) {

  // Examine the given dimensions
  if (width != free_number_format && abs(width) >= 62) {
    rtErr("The requested number exceeds format limits (maximum 64 characters, " +
          std::to_string(width) + " requested).", "intToString");
  }

  char buffer[64];
  if (width != free_number_format) {
    switch (style) {
    case NumberPrintStyle::STANDARD:
      snprintf(buffer, 64, "%*lld", width, value);
      break;
    case NumberPrintStyle::LEADING_ZEROS:
      snprintf(buffer, 64, "%0*lld", width, value);
      break;
    }
  }
  else {
    snprintf(buffer, 64, "%lld", value);
  }
  return std::string(buffer);
}

//-------------------------------------------------------------------------------------------------
const std::vector<TextGuard> operator+(const std::vector<TextGuard> &lhs,
                                       const std::vector<TextGuard> &rhs) {
  std::vector<TextGuard> result(lhs.begin(), lhs.end());
  result.insert(result.end(), rhs.begin(), rhs.end());
  return result;
}

//-------------------------------------------------------------------------------------------------
bool detectGuard(const TextFileReader &tfr, const int line_idx, const int pos_idx,
                 const std::string &guard_seq) {
  const int gsq_length = guard_seq.size();
  if (gsq_length == 0) {
    return false;
  }
  if (tfr.line_limits[line_idx] + pos_idx + gsq_length <= tfr.line_limits[line_idx + 1]) {
    bool matched = true;
    const int readbase = tfr.line_limits[line_idx] + pos_idx;
    for (int i = 0; i < gsq_length; i++) {
      matched = (matched && guard_seq[i] == tfr.text[readbase + i]);
    }
    return matched;
  }
  return false;
}

//-------------------------------------------------------------------------------------------------
bool detectGuard(const char* line, const int pos_idx, const std::string &guard_seq) {
  const int gsq_length = guard_seq.size();
  if (gsq_length == 0) {
    return false;
  }
  bool matched = true;
  for (int i = 0; i < gsq_length; i++) {
    matched = (matched && guard_seq[i] == line[pos_idx + i]);
  }
  return matched;
}

//-------------------------------------------------------------------------------------------------
bool detectGuard(const char* line, const size_t pos_idx, const size_t line_limit,
                 const std::string &guard_seq) {
  const size_t gsq_length = guard_seq.size();
  if (gsq_length == 0 || pos_idx + gsq_length > line_limit) {
    return false;
  }
  bool matched = true;
  for (size_t i = 0; i < gsq_length; i++) {
    matched = (matched && guard_seq[i] == line[pos_idx + i]);
  }
  return matched;
}

//-------------------------------------------------------------------------------------------------
int applyGuard(const char* line, const size_t pos_idx, const size_t line_limit,
               const std::vector<TextGuard> &markers) {
  const int n_markers = markers.size();
  for (int i = 0; i < n_markers; i++) {
    if (detectGuard(line, pos_idx, line_limit, markers[i].getLeft())) {
      return i;
    }
  }
  return -1;
}

//-------------------------------------------------------------------------------------------------
int applyGuard(const TextFileReader &tfr, const int line_idx, const int pos_idx,
               const std::vector<TextGuard> &markers) {
  const int n_markers = markers.size();
  for (int i = 0; i < n_markers; i++) {
    if (detectGuard(tfr, line_idx, pos_idx, markers[i].getLeft())) {
      return i;
    }
  }
  return -1;
}

//-------------------------------------------------------------------------------------------------
std::vector<bool> markEscapedCharacters(const char* textstr, const size_t n_char,
                                        const std::vector<TextGuard> &escapes) {
  const int n_escape = escapes.size();
  std::vector<bool> result(n_char, false);
  for (size_t i = 0; i < n_char; i++) {
    for (int j = 0; j < n_escape; j++) {
      if (detectGuard(textstr, i, escapes[j].getLeft())) {
        result[i + escapes[j].leftSize()] = true;
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<bool> markEscapedCharacters(const std::string &textstr,
                                        const std::vector<TextGuard> &escapes) {
  return markEscapedCharacters(textstr.data(), textstr.size(), escapes);
}

//-------------------------------------------------------------------------------------------------
std::vector<bool> markEscapedCharacters(const TextFileReader &tfr,
                                        const std::vector<TextGuard> &escapes) {
  std::vector<bool> result(tfr.line_limits[tfr.line_count], false);
  for (int i = 0; i < tfr.line_count; i++) {
    for (int j = tfr.line_limits[i]; j < tfr.line_limits[i + 1]; j++) {
      int guarded = applyGuard(tfr, i, j - tfr.line_limits[i], escapes);
      if (guarded >= 0 && j + escapes[guarded].leftSize() < tfr.line_limits[i + 1]) {
        result[j + escapes[guarded].leftSize()] = true;
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<bool> markGuardedText(const char* textstr, const size_t n_char,
                                  const size_t* line_limits, const int line_count,
                                  const std::vector<TextGuard> &markers,
                                  const std::vector<TextGuard> &alternatives,
                                  const std::vector<TextGuard> &escapes) {

  // Make all characters affected by an escape sequence.  Escapes are TextGuards with only a
  // left-hand sequence (by default, '\') and affect only the character directly behind them.
  std::vector<bool> literals = markEscapedCharacters(textstr, n_char, escapes);

  // Make an array of all markers containing the primary markers followed by the secondaries.
  // Record the number of primary markers, as these will take priority.
  std::vector<TextGuard> all_markers(markers);
  const int n_primary_markers = all_markers.size();
  all_markers.insert(all_markers.end(), alternatives.begin(), alternatives.end());
  int guarded = -1;
  std::vector<bool> result(n_char, false);
  for (int i = 0; i < line_count; i++) {
    for (size_t j = line_limits[i]; j < line_limits[i + 1]; j++) {

      // Process literals to have the same guarded state the escape sequences that preceded them
      if (literals[j]) {
        result[j] = result[j-1];
        continue;
      }

      // Check for any single-character text guard markers
      if (guarded >= 0) {
        if (detectGuard(textstr, j, line_limits[i + 1], all_markers[guarded].getRight())) {
          if (guarded < n_primary_markers) {

            // Characters in the guard sequence count as part of the guarded section of text,
            // and this section was guarded in one of the manners we want to detect.
            for (int k = 0; k < all_markers[guarded].rightSize(); k++) {
              result[j + k] = true;
            }
          }
          j += all_markers[guarded].rightSize() - 1;
          guarded = -1;
        }
        else {
          if (guarded < n_primary_markers) {
            result[j] = true;
          }
        }
      }
      else {

        // Try detecting a right guard in the absence of a guarded state.  If the left guard
        // of that guard pair is the same as the right guard, it doesn't count (false alarm),
        // as it means a sequence is just beginning.
        for (int k = 0; k < n_primary_markers; k++) {
          if (detectGuard(textstr, j, line_limits[i + 1], all_markers[k].getRight()) &&
              detectGuard(textstr, j, line_limits[i + 1], all_markers[k].getLeft()) == false) {
            rtErr("A mismatched right-hand text guard was detected on line " + std::to_string(i) +
                  ".", "markGuardedText");
          }
        }

        // Try detecting a guard
        guarded = applyGuard(textstr, j, line_limits[i + 1], all_markers);
        if (guarded >= 0) {

          // A guarded sequence of one of the types we want to detect just began.  The left
          // delimiter is a part of the sequence.
          if (guarded < n_primary_markers) {
            for (int k = 0; k < all_markers[guarded].leftSize(); k++) {
              result[j + k] = true;
            }
          }
          j += all_markers[guarded].leftSize();

          // If the guard extends to the end of the line, take it there and no further
          if (all_markers[guarded].rightSize() == 0) {
            if (guarded < n_primary_markers) {
              for (int k = j; k < line_limits[i + 1]; k++) {
                result[k] = true;
              }
            }
            j = line_limits[i + 1];
            guarded = -1;
          }

          // Decrement j before the loop iteration increments it
          j--;
        }
      }
    }

    // If we reach the end of a line and remain in a guarded state, make sure that the guard
    // can span multiple lines.
    if (guarded >= 0 && all_markers[guarded].getTerminationRequirement() &&
        all_markers[guarded].getSpan() != LineSpan::MULTIPLE) {
      rtErr("Guard sequence (\"" + all_markers[guarded].getLeft() + "\", \"" +
            all_markers[guarded].getRight() + "\") begins on line " + std::to_string(i + 1) +
            ".  It requires a terminating sequence and cannot span multiple lines.",
            "markGuardedText");
    }
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
std::vector<bool> markGuardedText(const char* textstr, const size_t n_char,
                                  const std::vector<TextGuard> &markers,
                                  const std::vector<TextGuard> &alternatives,
                                  const std::vector<TextGuard> &escapes) {

  // Create a catalog of points at which line breaks occur, after each carriage return.  Include
  // line limits at the front of the string, and at the back if the final character is not a
  // carriage return.
  std::vector<size_t> line_limits;
  size_t ncr = 1;
  for (size_t i = 0; i < n_char; i++) {
    if (textstr[i] == '\n') {
      ncr++;
    }
  }
  if (textstr[n_char - 1] != '\n') {
    ncr++;
  }
  line_limits.reserve(ncr);
  line_limits.push_back(0);
  for (size_t i = 0; i < n_char; i++) {
    if (textstr[i] == '\n') {
      line_limits.push_back(i + 1);
    }
  }
  if (textstr[n_char - 1] != '\n') {
    line_limits.push_back(n_char);
  }

  // Process the text with line demarcations using the same operations as a TextFile object.
  return markGuardedText(textstr, n_char, line_limits.data(), line_limits.size() - 1,
                         markers, alternatives, escapes);
}

//-------------------------------------------------------------------------------------------------
std::vector<bool> markGuardedText(const std::string &str,
                                  const std::vector<TextGuard> &markers,
                                  const std::vector<TextGuard> &alternatives,
                                  const std::vector<TextGuard> &escapes) {
  return markGuardedText(str.data(), str.size(), markers, alternatives, escapes);
}

//-------------------------------------------------------------------------------------------------
std::vector<bool> markGuardedText(const TextFileReader &tfr,
                                  const std::vector<TextGuard> &markers,
                                  const std::vector<TextGuard> &alternatives,
                                  const std::vector<TextGuard> &escapes) {
  return markGuardedText(tfr.text, tfr.line_limits[tfr.line_count], tfr.line_limits,
                         tfr.line_count, markers, alternatives, escapes);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> resolveScopes(const char* input_text, int length,
                               const std::vector<TextGuard> scope_guards,
                               const ExceptionResponse policy) {
  const int n_char = (length == 0) ? strlen(input_text) : length;
  int level = 0;
  std::vector<int> result(n_char, level);
  const int n_guard = scope_guards.size();
  std::vector<int> guard_layers;
  int i = 0;
  while (i < n_char) {

    // Check for any opening guard
    bool guard_found = false;
    for (int j = 0; j < n_guard; j++) {
      if (detectGuard(input_text, i, scope_guards[j].getLeft())) {
        guard_layers.push_back(j);
        const int n_left = scope_guards[j].leftSize();
        for (int k = 0; k < n_left; k++) {
          result[i + k] = level;
        }
        level++;
        i += n_left;
        guard_found = true;
        break;
      }
    }
    if (guard_found) {
      continue;
    }

    // Check for any closing guard 
    for (int j = 0; j < n_guard; j++) {
      if (detectGuard(input_text, i, scope_guards[j].getRight())) {
        if (j == guard_layers.back()) {
          level--;
          guard_layers.pop_back();
          const int n_right = scope_guards[j].rightSize();
          for (int k = 0; k < n_right; k++) {
            result[i + k] = level;
          }
          i += n_right;
          guard_found = true;
          break;
        }
        else {
          rtWarn("Closing guard \"" + scope_guards[j].getRight() + "\" found when expecting \"" +
                 scope_guards[guard_layers.back()].getRight() + "\" to close \"" +
                 scope_guards[guard_layers.back()].getLeft() + "\".", "resolveScopes");
        }
      }
    }
    if (guard_found) {
      continue;
    }

    // Final possibility: mark the one character as being at the current level and increment
    // the counter by one.
    result[i] = level;
    i++;
  }

  // Check that all scopes have been resolved
  if (guard_layers.size() > 0) {
    std::string all_remaining_scopes;
    const int n_rem = guard_layers.size();
    for (int i = 0; i < n_rem; i++) {
      all_remaining_scopes += "\"" + scope_guards[guard_layers[i]].getLeft() + "\"";
      if (i < n_rem - 1) {
        all_remaining_scopes += ", ";
      }
    }
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("String " + std::string(input_text) + " leaves " +
            std::to_string(guard_layers.size()) + " scopes unresolved: " + all_remaining_scopes +
            ".", "resolvedScopes");
    case ExceptionResponse::WARN:
      rtWarn("String " + std::string(input_text) + " leaves " +
             std::to_string(guard_layers.size()) + " scopes unresolved: " + all_remaining_scopes +
             ".", "resolvedScopes");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }

  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> resolveScopes(const std::string &input_text,
                               const std::vector<TextGuard> scope_guards,
                               const ExceptionResponse policy) {
  return resolveScopes(input_text.c_str(), input_text.size(), scope_guards, policy);
}

//-------------------------------------------------------------------------------------------------
int countDelimiters(const char* text, const size_t start_pos, const size_t length,
                    const char* delimiters, const int ndelim, const std::vector<bool> &quote_mask,
                    const std::vector<bool> &comment_mask) {
  int result = 0;
  const bool has_quote_mask   = (quote_mask.size() >= start_pos + length);
  const bool has_comment_mask = (comment_mask.size() >= start_pos + length);
  for (size_t i = 0; i < length; i++) {
    const size_t sp_i = start_pos + i;
    if ((has_quote_mask && quote_mask[sp_i]) || (has_comment_mask && comment_mask[sp_i])) {
      continue;
    }
    const char tmpc = text[sp_i];
    bool is_delimiter = false;
    for (int j = 0; j < ndelim; j++) {
      is_delimiter = (is_delimiter || (tmpc == delimiters[j]));
    }
    if (is_delimiter) {
      result++;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int countDelimiters(const std::string &text, const std::vector<char> &delimiters,
                    const std::vector<bool> &quote_mask, const std::vector<bool> &comment_mask) {
  return countDelimiters(text.data(), 0, text.size(), delimiters.data(), delimiters.size(),
                         quote_mask, comment_mask);
}

//-------------------------------------------------------------------------------------------------
int countWords(const char* line, const int start_pos, const int length) {
  const int end_pos = start_pos + length;
  int word_count = 0;
  bool on_word = false;
  for (int i = start_pos; i < end_pos; i++) {
    if (line[i] != ' ') {
      if (on_word == false) {
        on_word = true;
      }
    }
    else {
      if (on_word) {
        on_word = false;
        word_count++;
      }
    }
  }
  if (on_word) {
    word_count++;
  }
  return word_count;
}

//-------------------------------------------------------------------------------------------------
int countWords(const std::string &line, const int start_pos, const int length) {
  const int line_length = line.size();
  if (start_pos >= line_length || start_pos + length >= line_length) {
    return 0;
  }
  else {
    return countWords(line.c_str(), start_pos, length);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int countWords(const TextFile &tf, const int line, const int start_pos, const int length) {
  if (line >= tf.getLineCount()) {
    return 0;
  }
  const int line_length = tf.getLineLength(line);
  if (start_pos >= line_length || start_pos + length >= line_length) {
    return 0;
  }
  else {
    return countWords(tf.getLinePointer(line), start_pos, length);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void digestText(const char* text, const size_t start_pos, const size_t length,
                const char* separators, const int n_separators, const char* delimiters,
                const int n_delimiters, char* output, const std::vector<bool> &quote_mask,
                const std::vector<bool> &comment_mask, std::vector<bool> *output_quoted,
                std::vector<bool> *output_commented) {

  // Iterate over the output positions, asuming that the output begins at array position 0.
  size_t i = 0;
  const size_t hlim = start_pos + length;
  const bool has_quote_mask = (quote_mask.size() >= hlim);
  const bool has_comment_mask = (comment_mask.size() >= hlim);
  const bool has_quote_output = (output_quoted != nullptr && output_quoted->size() >= hlim);
  const bool has_comment_output = (output_commented != nullptr &&
                                   output_commented->size() >= hlim);
  for (size_t j = start_pos; j < hlim; j++) {
    if ((has_comment_mask && comment_mask[j]) || (has_quote_mask && quote_mask[j])) {
      
      // Because std::vector<bool> stores multiple boolean bits in one byte, it does not have a
      // .data() member variable and must be accessed using the at() member function.
      if (has_comment_output) {
        output_commented->at(i) = comment_mask[j];
      }
      if (has_quote_output) {
        output_quoted->at(i) = quote_mask[j];
      }
      output[i] = text[j];
      i++;
      continue;
    }
    bool is_delimiter = false;
    for (size_t k = 0; k < n_delimiters; k++) {
      is_delimiter = (is_delimiter || (text[j] == delimiters[k]));
    }
    if (is_delimiter) {
      output[i] = ' ';
      i++;
      continue;
    }
    bool is_separator = false;
    for (size_t k = 0; k < n_separators; k++) {
      is_separator = (is_separator || (text[j] == separators[k]));
    }
    if (is_separator) {
      output[i] = ' ';
      i++;
      output[i] = text[j];
      i++;
      output[i] = ' ';
      i++;
    }
    else {
      output[i] = text[j];
      i++;
    }
  }
}

//-------------------------------------------------------------------------------------------------
std::string digestText(const char* text, const size_t start_pos, const size_t length,
                       const char* separators, const int n_separators, const char* delimiters,
                       const int n_delimiters, const std::vector<bool> &quote_mask,
                       const std::vector<bool> &comment_mask, std::vector<bool> *output_quoted,
                       std::vector<bool> *output_commented) {
  
  // Count the number of separators.  The countDelimiters() function can be used in this respect.
  const int sep_instances = countDelimiters(text, start_pos, length, separators, n_separators,
                                            quote_mask, comment_mask);

  // Allocate the result to accommodate the number of separators identified.
  const size_t total_length = length + static_cast<size_t>(2 * sep_instances);
  std::string result(total_length, ' ');
  if (output_quoted != nullptr && output_quoted->size() < total_length) {
    output_quoted->resize(total_length);
  }
  if (output_commented != nullptr && output_commented->size() < total_length) {
    output_commented->resize(total_length);
  }
  digestText(text, start_pos, length, separators, n_separators, delimiters, n_delimiters,
             result.data(), quote_mask, comment_mask, output_quoted, output_commented);
  return result;
}

//-------------------------------------------------------------------------------------------------
void digestText(const std::string &text, const std::vector<char> &separators,
                const std::vector<char> &delimiters, std::string *output,
                const std::vector<bool> &quote_mask, const std::vector<bool> &comment_mask,
                std::vector<bool> *output_quoted, std::vector<bool> *output_commented) {

  // Count the number of separators in preparation for re-allocating the output if necessary.
  const int sep_instances = countDelimiters(text.data(), 0, text.size(), separators.data(),
                                            separators.size(), quote_mask, comment_mask);
  const size_t total_length = text.size() + static_cast<size_t>(2 * sep_instances);
  if (output->size() < total_length) {
    output->resize(total_length);
  }
  if (output_quoted != nullptr && output_quoted->size() < total_length) {
    output_quoted->resize(total_length);
  }
  if (output_commented != nullptr && output_commented->size() < total_length) {
    output_commented->resize(total_length);
  }
  digestText(text.data(), 0, text.size(), separators.data(), separators.size(), delimiters.data(),
             delimiters.size(), output->data(), quote_mask, comment_mask, output_quoted,
             output_commented);
}

//-------------------------------------------------------------------------------------------------
std::string digestText(const std::string &text, const std::vector<char> &separators,
                       const std::vector<char> &delimiters, const std::vector<bool> &quote_mask,
                       const std::vector<bool> &comment_mask, std::vector<bool> *output_quoted,
                       std::vector<bool> *output_commented) {
  const int sep_instances = countDelimiters(text.data(), 0, text.size(), separators.data(),
                                            separators.size(), quote_mask, comment_mask);
  const size_t total_length = text.size() + static_cast<size_t>(2 * sep_instances);
  std::string result(total_length, ' ');
  if (output_quoted != nullptr && output_quoted->size() < total_length) {
    output_quoted->resize(total_length);
  }
  if (output_commented != nullptr && output_commented->size() < total_length) {
    output_commented->resize(total_length);
  }
  digestText(text.data(), 0, text.size(), separators.data(), separators.size(), delimiters.data(),
             delimiters.size(), result.data(), quote_mask, comment_mask, output_quoted,
             output_commented);
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> separateText(const char* text, int n_char,
                                      const std::vector<bool> &comment_mask,
                                      const std::vector<bool> &quotation_mask,
                                      const std::vector<TextGuard> &quote_marks,
                                      const std::vector<std::string> &delimiters,
                                      const std::vector<TextGuard> &escapes) {

  // Check inputs
  const bool comments_exist = (comment_mask.size() > 0);
  const bool quotes_exist = (quotation_mask.size() > 0);
  const bool escapes_exist = (escapes.size() > 0);
  const int n_quote_guards = quote_marks.size();
  if (quotes_exist && n_quote_guards == 0) {
    rtErr("Separating text with masked quotations requires a list of quotation symbols.",
          "separateText");
  }
  if ((comments_exist && comment_mask.size() != n_char) ||
      (quotes_exist && quotation_mask.size() != n_char)) {
    std::string sample;
    for (int i = 0; i < std::min(n_char, 24); i++) {
      sample += text[i];
    }
    const std::string problem_area = (comments_exist && comment_mask.size() != n_char) ?
                                     "comment" : "quotation";
    rtErr("The " + problem_area + " mask for text beginning \"" + sample +
          "\" does not conform to the overall size of the text.", "separateText");
  }

  // Mark escaped characters
  const std::vector<bool> escape_mask = (escapes_exist) ?
                                        markEscapedCharacters(text, n_char, escapes) :
                                        std::vector<bool>(n_char, false);

  // Convert all additional delimiters to whitespace
  const int n_delimiters = delimiters.size();
  std::vector<bool> delimiter_mask(n_char, false);
  for (int i = 0; i < n_char; i++) {
    delimiter_mask[i] = (text[i] == ' ');
  }
  for (int i = 0; i < n_delimiters; i++) {
    const int delimiter_size = delimiters[i].size();
    if (delimiter_size == 1) {
      const char idel = delimiters[i][0];
      if (comments_exist && quotes_exist) {
        for (int j = 0; j < n_char; j++) {
          delimiter_mask[j] = ((comment_mask[j] == false && quotation_mask[j] == false &&
                                escape_mask[j] == false && text[j] == idel) || delimiter_mask[j]);
        }
      }
      else if (comments_exist) {
        for (int j = 0; j < n_char; j++) {
          delimiter_mask[j] = ((comment_mask[j] == false && escape_mask[j] == false &&
                                text[j] == idel) || delimiter_mask[j]);
        }
      }
      else if (quotes_exist) {
        for (int j = 0; j < n_char; j++) {
          delimiter_mask[j] = ((quotation_mask[j] == false && escape_mask[j] == false &&
                                text[j] == idel) || delimiter_mask[j]);
        }
      }
      else {
        for (int j = 0; j < n_char; j++) {
          delimiter_mask[j] = ((escape_mask[j] == false && text[j] == idel) || delimiter_mask[j]);
        }
      }
    }
    else {
      const char* idel = delimiters[i].c_str();
      for (int j = 0; j < n_char; j++) {
        const int dlim = j + delimiters[i].size();
        if (dlim >= n_char) {
          continue;
        }
        const char* idel = delimiters[i].c_str();
        bool mismatch = false;
        if (comments_exist && quotes_exist) {
          for (int k = j; k < dlim; k++) {
            mismatch = ((comment_mask[k] == false && quotation_mask[k] == false &&
                         escape_mask[k] == false && text[k] != idel[k - j]) || mismatch);
          }
        }
        else if (comments_exist) {
          for (int k = j; k < dlim; k++) {
            mismatch = ((comment_mask[k] == false && escape_mask[k] == false &&
                         text[k] != idel[k - j]) || mismatch);
          }
        }
        else if (quotes_exist) {
          for (int k = j; k < dlim; k++) {
            mismatch = ((quotation_mask[k] == false && escape_mask[k] == false &&
                         text[k] != idel[k - j]) || mismatch);
          }
        }
        else {
          for (int k = j; k < dlim; k++) {
            mismatch = ((escape_mask[k] == false && text[k] != idel[k - j]) || mismatch);
          }
        }
        if (mismatch == false) {
          for (int k = j; k < dlim; k++) {
            delimiter_mask[k] = true;
          }
        }
      }
    }
  }

  // With the delimiter mask, comments mask, quotations mask, and
  // escaped characters mask in hand, find distinct strings.
  std::string growing_word;
  std::vector<std::string> result;
  for (int i = 0; i < n_char; i++) {

    // As escaped character just gets added to the growing word string, preventing anything
    // else from happening.
    if (escape_mask[i]) {
      growing_word += text[i];
      continue;
    }

    // Continue until a non-delimiter, non-comment character is encountered
    if (delimiter_mask[i] || (comments_exist && comment_mask[i])) {
      if (growing_word.size() > 0) {
        result.push_back(growing_word);
        growing_word.resize(0);
      }
      continue;
    }

    // A quotation is one word, unless it's multiple quotes back-to-back
    if (quotes_exist) {
      int n_quoted_char = 0;
      while (i + n_quoted_char < n_char && quotation_mask[i + n_quoted_char]) {
        n_quoted_char++;
      }
      if (n_quoted_char > 0) {

        // Search all delimiters and determine which quotes are responsible for this string
        int guarded = -1;
        for (int j = i; j < i + n_quoted_char; j++) {

          // Again, escaped characters just get added and pass over any further processing.
          // The guarded state (some sort of quotation) does not change.
          if (escape_mask[j] && guarded >= 0) {
            growing_word += text[j];
            continue;
          }

          if (guarded == -1) {
            for (int k = 0; k < n_quote_guards; k++) {
              if (j + quote_marks[k].leftSize() < i + n_quoted_char &&
                  detectGuard(text, j, quote_marks[k].getLeft())) {
                guarded = k;
                j += quote_marks[k].leftSize();
                break;
              }
            }
          }

          // Some form of guard must be in place by now
          if (guarded < 0) {
            std::string sample;
            for (int k = 0; k < std::min(n_quoted_char, 240); k++) {
              sample += text[i + k];
            }
            rtErr("Quotation beginning <<" + sample + ">> does not conform to any of the quote "
                  "marks provided.", "separateText");
          }

          // Check for a right-hand quotation terminator
          if (j + quote_marks[guarded].rightSize() <= i + n_quoted_char &&
              detectGuard(text, j, quote_marks[guarded].getRight())) {
            if (growing_word.size() > 0) {
              result.push_back(growing_word);
              growing_word.resize(0);
            }
            j += quote_marks[guarded].rightSize() - 1;
            guarded = -1;
          }

          // Grow the current word if this is still under a quotation guard
          if (guarded >= 0) {
            growing_word += text[j];
          }
        }

        // The quotation should end with a closing mark
        if (guarded >= 0) {
          std::string sample;
          for (int k = 0; k < std::min(n_quoted_char, 24); k++) {
            sample += text[i + k];
          }
          rtErr("Quotation beginning <<" + sample + ">> terminates without using any of the quote "
                "marks provided.", "separateText");
        }

        // Increment the counter to reflect progress through this quoted stretch of text
        i += n_quoted_char - 1;
      }
    }
    if (delimiter_mask[i] == false && (comments_exist == false || comment_mask[i] == false) &&
        (quotes_exist == false || quotation_mask[i] == false)) {

      // Grow the current word
      growing_word += text[i];
    }
  }

  // Append the final word 
  if (growing_word.size() > 0) {
    result.push_back(growing_word);
  }

  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> separateText(const std::string &text,
                                      const std::vector<bool> &comment_mask,
                                      const std::vector<bool> &quotation_mask,
                                      const std::vector<TextGuard> &quote_marks,
                                      const std::vector<std::string> &delimiters,
                                      const std::vector<TextGuard> &escapes) {
  return separateText(text.c_str(), text.size(), comment_mask, quotation_mask, quote_marks,
                      delimiters, escapes);
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> separateText(const std::string &text,
                                      const std::vector<std::string> &delimiters,
                                      const std::vector<TextGuard> &escapes) {
  return separateText(text.c_str(), text.size(), std::vector<bool>(), std::vector<bool>(),
                      std::vector<TextGuard>(), delimiters, escapes);
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> separateText(const TextFile &tf,
                                      const std::vector<bool> &comment_mask,
                                      const std::vector<bool> &quotation_mask,
                                      const std::vector<TextGuard> &quote_marks,
                                      const std::vector<std::string> &delimiters,
                                      const std::vector<TextGuard> &escapes) {
  TextFileReader tfr = tf.data();
  return separateText(tfr.text, tfr.line_limits[tfr.line_count], comment_mask, quotation_mask,
                      quote_marks, delimiters, escapes);
}
  
//-------------------------------------------------------------------------------------------------
int findStringInVector(const std::vector<std::string> &vec, const std::string &query) {

  // Check for simple problems
  const size_t qsize = query.size();
  const int vsize = vec.size();
  if (vsize == 0) {
    return 0;
  }

  // Peek at the first few characters to see if a match is possible.
  // There is at least one character in the query string.  Do a deeper search if warranted.
  const char first_char = query[0];  
  for (int i = 0; i < vsize; i++) {
    if (qsize == 0LLU && vec[i].size() == qsize) {
      return i;
    }
    if (vec[i].size() == qsize && vec[i][0] == first_char && vec[i] == query) {
      return i;
    }
  }

  // The query was not found.  Return the length of the vector.
  return vsize;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> vectorStrtol(const std::vector<std::string> &sv, const ExceptionResponse policy) {
  const int nval = sv.size();
  std::vector<int> result(nval, 0.0);
  for (int i = 0; i < nval; i++) {
    if (verifyNumberFormat(sv[i].c_str(), NumberFormat::INTEGER)) {
      result[i] = strtol(sv[i].c_str(), nullptr, 10);
    }
    else {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("The string " + sv[i] + " is not a valid integer.", "vectorStrtol");
      case ExceptionResponse::WARN:
        rtWarn("The string " + sv[i] + " is not a valid integer.  Substituting 0.",
               "vectorStrotol");
	break;
      case ExceptionResponse::SILENT:
	break;
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> vectorStrtod(const std::vector<std::string> &sv,
                                 const ExceptionResponse policy) {
  const int nval = sv.size();
  std::vector<double> result(nval, 0.0);
  for (int i = 0; i < nval; i++) {
    if (verifyNumberFormat(sv[i].c_str(), NumberFormat::STANDARD_REAL) ||
        verifyNumberFormat(sv[i].c_str(), NumberFormat::SCIENTIFIC)) {
      result[i] = strtod(sv[i].c_str(), nullptr);
    }
    else {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("The string " + sv[i] + " is not a valid real number.", "vectorStrtod");
      case ExceptionResponse::WARN:
        rtWarn("The string " + sv[i] + " is not a valid real number.  Substituting 0.0.",
               "vectorStrtod");
	break;
      case ExceptionResponse::SILENT:
	break;
      }
    }
  }
  return result;
}

} // namespace parse
} // namespace stormm
