#include <climits>
#include "copyright.h"
#include "Reporting/error_format.h"
#include "ascii_numbers.h"

namespace stormm {
namespace parse {

//-------------------------------------------------------------------------------------------------
void printNumberSeries(std::ofstream *foutp, const std::vector<PolyNumeric> &values,
                       const int values_per_line, const int width, const int decimal,
                       const NumberFormat format, const std::string &caller,
                       const std::string &task) {
  
  // Allocate an individual line
  const int nl_char = (values_per_line * width) + 65;
  std::vector<char> line(nl_char, '\0');

  // Check integer format for expedited ASCII printing
  const bool fixed_resources_exceeded = (sizeof(int) != 4 || decimal > 9);

  // Codify the numbers and symbols that might be printed.
  const int e_index = 'e';
  const int zero_index = '0';
  const int space_index = ' ';
  const int positive_index = '+';
  const int negative_index = '-';
  const int nz_index = negative_index - zero_index;
  const int ns_index = negative_index - space_index;
  const int pn_index = positive_index - negative_index;
  const int zs_index = zero_index - space_index;
  const int dot_index = '.';

  // Make a vector of integers to prevent type recasting until the line is completely ready.
  std::vector<int> line_symbols(width * values_per_line, space_index);
  std::vector<int> sign_symbols(values_per_line, space_index);
  std::vector<double> tmp_values(values_per_line, 0.0);
  std::vector<int> integral(values_per_line, 0);
  std::vector<int> integral_nonzero(values_per_line, 0);
  std::vector<int> last_blank_char(values_per_line, 0);
  std::vector<int> fraction(values_per_line, 0);

  // Some initialization is common to all switch cases
  const int nval = values.size();
  int llim = 0;

  // Real numbers and integers can all be printed as number series
  switch (format) {
  case NumberFormat::SCIENTIFIC:
    {
      const int dec_offset = width - decimal - 5;
      const int exp_offset = width - 4;
      const int sign_offset = width - decimal - 7;

      // Loop over all values.  Filling out scientific notation is easier than general format.
      while (llim < nval) {

        // Nested loops proceed over all values on the present line
        const int hlim = std::min(llim + values_per_line, nval);
        const int n_this_line = hlim - llim;

        // Check for extreme values that would break the format
        bool broken_values = false;
        for (int j = llim; j < hlim; j++) {
          const int base_ten_exp = log10(fabs(values[j].d));
          broken_values = (base_ten_exp < -99 || base_ten_exp > 99 || broken_values);
        }
        if (broken_values) {
          for (int j = llim; j < hlim; j++) {
            const int base_ten_exp = (fabs(values[j].d) > 1.0e-98) ? log10(fabs(values[j].d)) : 0;
            if (base_ten_exp < -99 || base_ten_exp > 99) {
              rtErr("A value of " + realToString(values[j].d, width, decimal, format) +
                    " cannot be represented in format %" + std::to_string(width) + "." +
                    std::to_string(decimal) + "e. (" + task + ".)", "printNumberSeries");
            }
          }
        }

        // Print the numbers
        for (int j = 0; j < n_this_line; j++) {
          snprintf(&line[j * width], nl_char - (j * width), "%*.*e", width, decimal,
                   values[llim + j].d);
        }
        line[n_this_line * width] = '\n';
        line[n_this_line * width + 1] = '\0';
        foutp->write(line.data(), n_this_line * width + 1);

        // Increment the lower value limit
        llim += values_per_line;        
      }
    }
    break;
  case NumberFormat::STANDARD_REAL:
    {
      // Assess the limits of the ascii text format
      const double pos_format_limit = pow(10.0, width - decimal - 1);
      const double neg_format_limit = pos_format_limit / -10.0;
    
      // Placement and sizing pre-computations for general format real numbers
      const double decimal_scale = pow(10.0, decimal);
      const int dec_offset = width - decimal - 1;
      if (width - decimal <= 1) {
        rtErr("Real numbers cannot be correctly printed with " + std::to_string(decimal) +
              " digits after the decimal and " + std::to_string(width) + " total format space. (" +
              task + ".)", "printNumberSeries");
      }
      const int integral_base_factor = static_cast<int>(pow(10.0, dec_offset - 1) + 0.01);
      const int fraction_base_factor = static_cast<int>(pow(10.0, decimal - 1) + 0.01);

      // Loop over all values
      while (llim < nval) {

        // Nested loops proceed over all values on the present line
        const int hlim = std::min(llim + values_per_line, nval);
        const int n_this_line = hlim - llim;

        // Check against numbers that are just too big for a 32-bit integer format, and trap
        // if numbers exceed the printed column format.
        bool extreme_values = fixed_resources_exceeded;
        bool broken_values = false;
        for (int j = llim; j < hlim; j++) {
          broken_values = (values[j].d >= pos_format_limit || values[j].d <= neg_format_limit ||
                           broken_values);
          extreme_values = (values[j].d > 2147483647.0 || values[j].d < -2147483647.0 ||
                            extreme_values);
        }
        if (broken_values) {
          for (int j = llim; j < hlim; j++) {
            if (values[j].d >= pos_format_limit || values[j].d <= neg_format_limit) {
              rtErr("A value of " + realToString(values[j].d, decimal) + " cannot be represented "
                    "in format %" + std::to_string(width) + "." + std::to_string(decimal) + "lf. "
                    "(" + task + ".)", "printNumberSeries");
            }
          }
        }

        // If the values are extreme (breaking the integer format before or after the decimal)
        // but not so much that they break the actual format, print them with the standard tools
        if (extreme_values) {
          for (int j = llim; j < hlim; j++) {
            snprintf(&line[(j - llim) * width], nl_char - ((j - llim) * width), "%*.*lf", width,
                     decimal, values[j].d);
          }
          line[n_this_line * width] = '\n';
          line[n_this_line * width + 1] = '\0';
          foutp->write(line.data(), n_this_line * width + 1);

          // Increment the lower value limit and continue
          llim += values_per_line;
          continue;
        }

        // If the values are all amenable to the format, print them with an expedited technique.
        for (int j = 0; j < n_this_line * width; j++) {
          line_symbols[j] = space_index;
        }
        for (int j = 0; j < n_this_line; j++) {
          tmp_values[j] = values[j + llim].d;
          sign_symbols[j] = space_index + (ns_index * (tmp_values[j] < 0.0));
          tmp_values[j] = fabs(tmp_values[j]);
          integral[j] = static_cast<int>(floor(tmp_values[j]) + 0.01);
          fraction[j] = static_cast<int>(round((tmp_values[j] -
                                                floor(tmp_values[j])) * decimal_scale) + 0.01);
          if (fraction[j] == decimal_scale) {
            fraction[j] = 0;
            integral[j] += 1;
          }
          line_symbols[j*width + dec_offset] = dot_index;
          integral_nonzero[j] = (integral[j] != 0);
          last_blank_char[j] = -1;
        }
        int kcon = 0;
        for (int k = integral_base_factor; k > 0; k /= 10) {
          for (int j = 0; j < n_this_line; j++) {
            const int integral_over_k = integral[j] / k;
            last_blank_char[j] += (integral_over_k == 0 && last_blank_char[j] == kcon - 1);
            line_symbols[j*width + kcon] += (last_blank_char[j] != kcon) * integral_nonzero[j] *
                                            (zs_index + integral_over_k);
            integral[j] -= integral_over_k * k;
          }
          kcon++;
        }
        kcon = dec_offset + 1;
        for (int k = fraction_base_factor; k > 0; k /= 10) {
          for (int j = 0; j < n_this_line; j++) {
            const int fraction_over_k = fraction[j] / k;
            line_symbols[j*width + kcon] = zero_index + fraction_over_k;
            fraction[j] -= fraction_over_k * k;
          }
          kcon++;
        }
        for (int j = 0; j < n_this_line; j++) {
          if (integral_nonzero[j] == 0) {
            if (dec_offset > 1) {
              line_symbols[j*width + dec_offset - 1] = zero_index;
              line_symbols[j*width + dec_offset - 2] = sign_symbols[j];
            }
            else {

              // The remaining option is that dec_offset == 1
              line[j * width] = zero_index + (sign_symbols[j] == negative_index) * nz_index;
            }
          }
          else if (last_blank_char[j] >= 0) {
            line_symbols[j*width + last_blank_char[j]] = sign_symbols[j];
          }
        }
        for (int j = 0; j < n_this_line * width; j++) {
          line[j] = line_symbols[j];
        }
        line[n_this_line * width] = '\n';
        line[n_this_line * width + 1] = '\0';
        foutp->write(line.data(), n_this_line * width + 1);

        // Increment the lower value limit
        llim += values_per_line;
      }
    }
    break;
  case NumberFormat::INTEGER:
    {
      // Compute the integer limit of the values for the format
      int acc = 1;
      const int max_powers_of_ten = std::min(width - 1, 9);
      for (int j = 0; j < max_powers_of_ten; j++) {
        acc *= 10;
      }
      const int pos_format_limit = (width >= 10) ? INT_MAX : acc - 1;
      const int neg_format_limit = (width >= 11) ? INT_MIN : -acc / 10 + 1;

      // Loop over all values
      while (llim < nval) {

        // Nested loops proceed over all values on the present line
        const int hlim = std::min(llim + values_per_line, nval);
        const int n_this_line = hlim - llim;

        // Check for extreme values that would break the format
        bool broken_values = false;
        for (int j = llim; j < hlim; j++) {
          broken_values = (values[j].i > pos_format_limit || values[j].i < neg_format_limit ||
                           broken_values);
        }
        if (broken_values) {
          for (int j = llim; j < hlim; j++) {
            if (values[j].i > pos_format_limit || values[j].i < neg_format_limit) {
              rtErr("A value of " + std::to_string(values[j].i) + " cannot be represented in "
                    "format %" + std::to_string(width) + "d. (" + task + ".)",
                    "printNumberSeries");
            }
          }
        }
        for (int j = 0; j < n_this_line; j++) {
          snprintf(&line[j * width], nl_char - (j * width), "%*d", width, values[llim + j].i);
        }
        line[n_this_line * width] = '\n';
        line[n_this_line * width + 1] = '\0';
        foutp->write(line.data(), n_this_line * width + 1);

        // Increment the lower value limit
        llim += values_per_line;
      }
    }
    break;
  case NumberFormat::LONG_LONG_INTEGER:
    {
      // Compute the integer limit of the values for the format
      long long int acc = 1;
      const int max_powers_of_ten = std::min(width - 1, 18);
      for (int j = 0; j < max_powers_of_ten; j++) {
        acc *= 10LL;
      }
      const long long int pos_format_limit = (width >= 19) ? LLONG_MAX : acc - 1LL;
      const long long int neg_format_limit = (width >= 20) ? LLONG_MIN : -acc / 10LL + 1;

      // Loop over all values
      while (llim < nval) {

        // Nested loops proceed over all values on the present line
        const int hlim = std::min(llim + values_per_line, nval);
        const int n_this_line = hlim - llim;

        // Check for extreme values that would break the format
        bool broken_values = false;
        for (int j = llim; j < hlim; j++) {
          broken_values = (values[j].lli > pos_format_limit || values[j].lli < neg_format_limit ||
                           broken_values);
        }
        if (broken_values) {
          for (int j = llim; j < hlim; j++) {
            if (values[j].lli > pos_format_limit || values[j].lli < neg_format_limit) {
              rtErr("A value of " + std::to_string(values[j].lli) + " cannot be represented in "
                    "format %" + std::to_string(width) + "lld. (" + task + ".)",
                    "printNumberSeries");
            }
          }
        }
        for (int j = 0; j < n_this_line; j++) {
          snprintf(&line[j * width], nl_char - (j * width), "%*lld", width, values[llim + j].lli);
        }
        line[n_this_line * width] = '\n';
        line[n_this_line * width + 1] = '\0';
        foutp->write(line.data(), n_this_line * width + 1);

        // Increment the lower value limit
        llim += values_per_line;
      }
    }
    break;
  case NumberFormat::CHAR4:
    {
      // Loop over all values.  This is simple.
      while (llim < nval) {

        // Nested loops proceed over all values on the present line
        const int hlim = std::min(llim + values_per_line, nval);
        const int n_this_line = hlim - llim;
        for (int j = 0; j < n_this_line; j++) {
          snprintf(&line[j * 5], nl_char - (j * 5), "%c%c%c%c ", values[llim + j].c4.x,
                   values[llim + j].c4.y, values[llim + j].c4.z, values[llim + j].c4.w);
        }
        line[n_this_line * width] = '\n';
        line[n_this_line * width + 1] = '\0';
        foutp->write(line.data(), n_this_line * width + 1);

        // Increment the lower value limit
        llim += values_per_line;
      }
    }
    break;
  case NumberFormat::UNSIGNED_INTEGER:
    {
      // Compute the integer limit of the values for the format
      unsigned int acc = 1;
      const int max_powers_of_ten = std::min(width - 1, 9);
      for (int j = 0; j < max_powers_of_ten; j++) {
        acc *= 10U;
      }
      const unsigned int pos_format_limit = (width >= 10) ? UINT_MAX : acc - 1;
      
      // Loop over all values.  Filling out scientific notation is easier than general format.
      while (llim < nval) {

        // Nested loops proceed over all values on the present line
        const int hlim = std::min(llim + values_per_line, nval);
        const int n_this_line = hlim - llim;

        // Check for extreme values that would break the format
        bool broken_values = false;
        for (int j = llim; j < hlim; j++) {
          broken_values = (values[j].ui > pos_format_limit || broken_values);
        }
        if (broken_values) {
          for (int j = llim; j < hlim; j++) {
            if (values[j].ui > pos_format_limit) {
              rtErr("A value of " + std::to_string(values[j].ui) + " cannot be represented in "
                    "format %" + std::to_string(width) + "u. (" + task + ".)",
                    "printNumberSeries");
            }
          }
        }
        for (int j = 0; j < n_this_line; j++) {
          snprintf(&line[j * width], nl_char - (j * width), "%*u", width, values[llim + j].ui);
        }
        line[n_this_line * width] = '\n';
        line[n_this_line * width + 1] = '\0';
        foutp->write(line.data(), n_this_line * width + 1);

        // Increment the lower value limit
        llim += values_per_line;
      }
    }
    break;
  case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
    {
      // Compute the integer limit of the values for the format
      unsigned int acc = 1;
      const int max_powers_of_ten = std::min(width - 1, 9);
      for (int j = 0; j < max_powers_of_ten; j++) {
        acc *= 10U;
      }
      const unsigned int pos_format_limit = (width >= 10) ? UINT_MAX : acc - 1;

      // Loop over all values.  Filling out scientific notation is easier than general format.
      while (llim < nval) {

        // Nested loops proceed over all values on the present line
        const int hlim = std::min(llim + values_per_line, nval);
        const int n_this_line = hlim - llim;

        // Check for extreme values that would break the format
        bool broken_values = false;
        for (int j = llim; j < hlim; j++) {
          broken_values = (values[j].ulli > pos_format_limit || broken_values);
        }
        if (broken_values) {
          for (int j = llim; j < hlim; j++) {
            if (values[j].ui > pos_format_limit) {
              rtErr("A value of " + std::to_string(values[j].ulli) + " cannot be represented in "
                    "format %" + std::to_string(width) + "llu. (" + task + ".)",
                    "printNumberSeries");
            }
          }
        }
        for (int j = 0; j < n_this_line; j++) {
          snprintf(&line[j * width], nl_char - (j * width), "%*llu", width, values[llim + j].ulli);
        }
        line[n_this_line * width] = '\n';
        line[n_this_line * width + 1] = '\0';
        foutp->write(line.data(), n_this_line * width + 1);

        // Increment the lower value limit
        llim += values_per_line;
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void printNumberSeries(std::ofstream *foutp, const std::vector<PolyNumeric> &values,
                       const int values_per_line, const int width, const NumberFormat format,
                       const std::string &caller, const std::string &task) {
  printNumberSeries(foutp, values, values_per_line, width, 0, format, caller, task);
}

//-------------------------------------------------------------------------------------------------
std::vector<PolyNumeric> readNumberSeries(const TextFile &tf, const int start_line,
                                          const int n_values, const int values_per_line,
                                          const int width, const int decimal,
                                          const NumberFormat format, const std::string &caller,
                                          const std::string &task) {

  // Codify the numbers and symbols that might be printed.
  const int e_index = 'e';
  const int ecap_index = 'E';
  const int zero_index = '0';
  const int nine_index = '9';
  const int space_index = ' ';
  const int positive_index = '+';
  const int negative_index = '-';
  const int nz_index = negative_index - zero_index;
  const int ns_index = negative_index - space_index;
  const int pn_index = positive_index - negative_index;
  const int zs_index = zero_index - space_index;
  const int dot_index = '.';

  // Prepare to store the result
  const TextFileReader tfr = tf.data();
  std::vector<PolyNumeric> result;
  result.resize(n_values);
  std::vector<int> line_integral(values_per_line, 0);
  std::vector<ulint> line_uintegral(values_per_line, 0ul);
  std::vector<llint> line_llintegral(values_per_line, 0LL);
  std::vector<ullint> line_ullintegral(values_per_line, 0ull);
  std::vector<int> line_decimal(values_per_line, 0);
  std::vector<int> line_negative(values_per_line, 0);
  std::vector<int> line_exponent(values_per_line, 0);
  const double lowest_significant_digit = pow(0.1, decimal);

  // Pre-flight check: does this data conform to a format that can be speed-read?
  const int last_line_may_hold = tfr.line_limits[tfr.line_count] -
                                 tfr.line_limits[tfr.line_count - 1];
  if ((tfr.line_count - start_line - 1) * values_per_line + last_line_may_hold < n_values) {
    const std::string callerline = (caller.size() > 0) ? "  Called by: " + caller + "." : "";
    const std::string taskline = (task.size() > 0) ? "  Task: " + task + "." : "";
    rtErr("Too many values were requested starting at line " + std::to_string(start_line + 1) +
          ", at the rate of " + std::to_string(values_per_line) + " per line.  File " + 
          tf.getFileName() + " will end prematurely." + callerline + taskline, "readNumberSeries");
  }
  int iline = start_line;
  int llim = 0;
  bool format_understood = true;
  while (llim < n_values && format_understood) {
    const int n_this_line = (llim + values_per_line <= n_values) ? values_per_line :
                                                                   n_values - llim;
    const int readbase = tfr.line_limits[iline];
    switch (format) {
    case NumberFormat::SCIENTIFIC:
      for (int i = 0; i < n_this_line; i++) {
        const int number_end = readbase + ((i + 1) * width);
        const int tmpc = tfr.text[number_end - 1];
        format_understood = (format_understood && tmpc >= zero_index && tmpc <= nine_index);
        const int tmp2c = tfr.text[number_end - 3];
        format_understood = (format_understood && (tmp2c == positive_index ||
                                                   tmp2c == negative_index));
        const int tmp3c = tfr.text[number_end - 4];
        format_understood = (format_understood && (tmp3c == e_index || tmp3c == ecap_index));
        const int tmp4c = tfr.text[number_end - 5 - decimal];
        format_understood = (format_understood && (tmp4c == dot_index));
      }
      break;
    case NumberFormat::STANDARD_REAL:
      for (int i = 0; i < n_this_line; i++) {
        const int number_end = readbase + ((i + 1) * width);
        const int tmpc = tfr.text[number_end - 1];
        format_understood = (format_understood && tmpc >= zero_index && tmpc <= nine_index);
        const int tmp2c = tfr.text[number_end - decimal - 1];
        format_understood = (format_understood && (tmp2c == dot_index));
      }
      break;
    case NumberFormat::INTEGER:
    case NumberFormat::LONG_LONG_INTEGER:
    case NumberFormat::UNSIGNED_INTEGER:
    case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
      for (int i = 0; i < n_this_line; i++) {
        const int tmpc = tfr.text[readbase + ((i + 1) * width) - 1];
        format_understood = (format_understood && tmpc >= zero_index && tmpc <= nine_index);
      }
      break;
    case NumberFormat::CHAR4:
      break;
    }
    iline++;
    llim += n_this_line;
  }

  // Whether the format seems to be understood or not, check that each number is interpretable.
  // The fixed resources are unnecessary with the CHAR4 number format, so they are, by
  // construction, not exceeded.
  const bool fixed_resources_exceeded = (format != NumberFormat::CHAR4 &&
                                         (sizeof(int) != 4 || decimal > 9 ||
                                          (format == NumberFormat::STANDARD_REAL &&
                                           width - decimal - 1 > 9)));
  std::string buffer;
  if (fixed_resources_exceeded) {
    buffer.resize(width);
  }
  int max_signs;
  iline = start_line;
  llim = 0;
  switch (format) {
  case NumberFormat::SCIENTIFIC:
    max_signs = 2;
    break;
  case NumberFormat::STANDARD_REAL:
  case NumberFormat::INTEGER:
  case NumberFormat::LONG_LONG_INTEGER:
  case NumberFormat::UNSIGNED_INTEGER:
  case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
  case NumberFormat::CHAR4:
    max_signs = 1;
    break;
  }
  while (llim < n_values) {
    const int n_this_line = (llim + values_per_line <= n_values) ? values_per_line :
                                                                   n_values - llim;
    const int readbase = tfr.line_limits[iline];
    bool problem = false;
    switch (format) {
    case NumberFormat::SCIENTIFIC:
      for (int i = 0; i < n_this_line; i++) {
        int nsigns = 0;
        int nexps = 0;
        bool number_begins = false;
	int prev_tmpc = 0;
        for (int j = 0; j < width; j++) {
          const int tmpc = tfr.text[readbase + (i * width) + j];
          number_begins = (number_begins || tmpc != space_index);
          nsigns += (tmpc == negative_index || tmpc == positive_index);
          nexps += (tmpc == e_index || tmpc == ecap_index);
          problem = (problem || ((tmpc < zero_index || tmpc > nine_index) && tmpc != dot_index &&
                                 tmpc != space_index && tmpc != e_index && tmpc != ecap_index &&
                                 tmpc != negative_index && tmpc != positive_index) ||
                     nsigns > 2 || nexps > 1 || (number_begins && tmpc == space_index) ||
                     ((tmpc == negative_index || tmpc == positive_index) && j > 0 &&
                      prev_tmpc != space_index && prev_tmpc != e_index &&
                      prev_tmpc != ecap_index) ||
                     (tmpc == dot_index && j > 0 && prev_tmpc != space_index &&
                      (prev_tmpc < zero_index || prev_tmpc > nine_index)) ||
		     ((tmpc == e_index || tmpc == ecap_index) && (prev_tmpc < zero_index ||
                                                                  prev_tmpc > nine_index)));
	  prev_tmpc = tmpc;
        }
      }
      break;
    case NumberFormat::STANDARD_REAL:
      for (int i = 0; i < n_this_line; i++) {
        int nsigns = 0;
        bool number_begins = false;
	int prev_tmpc = 0;
        for (int j = 0; j < width; j++) {
          const int tmpc = tfr.text[readbase + (i * width) + j];
          number_begins = (number_begins || tmpc != space_index);
          nsigns += (tmpc == negative_index || tmpc == positive_index);
          problem = (problem || ((tmpc < zero_index || tmpc > nine_index) && tmpc != dot_index &&
                                 tmpc != space_index && tmpc != negative_index &&
                                 tmpc != positive_index) || nsigns > 1 ||
                     (number_begins && tmpc == space_index) ||
                     ((tmpc == negative_index || tmpc == positive_index) && j > 0 &&
                      prev_tmpc != space_index));
	  prev_tmpc = tmpc;
        }
      }
      break;
    case NumberFormat::INTEGER:
    case NumberFormat::LONG_LONG_INTEGER:
      for (int i = 0; i < n_this_line; i++) {
        int nsigns = 0;
        bool number_begins = false;
	int prev_tmpc = 0;
        for (int j = 0; j < width; j++) {
          const int tmpc = tfr.text[readbase + (i * width) + j];
          number_begins = (number_begins || tmpc != space_index);
          nsigns += (tmpc == negative_index || tmpc == positive_index);
          problem = (problem || ((tmpc < zero_index || tmpc > nine_index) && tmpc != space_index &&
                                 tmpc != negative_index && tmpc != positive_index) || nsigns > 1 ||
                     (number_begins && tmpc == space_index) ||
                     ((tmpc == negative_index || tmpc == positive_index) && j > 0 &&
                      prev_tmpc != space_index));
	  prev_tmpc = tmpc;
        }
      }
      break;
    case NumberFormat::UNSIGNED_INTEGER:
    case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
      for (int i = 0; i < n_this_line; i++) {
        int nsigns = 0;
        bool number_begins = false;
	int prev_tmpc = 0;
        for (int j = 0; j < width; j++) {
          const int tmpc = tfr.text[readbase + (i * width) + j];
          number_begins = (number_begins || tmpc != space_index);
          nsigns += (tmpc == positive_index);
          problem = (problem || ((tmpc < zero_index || tmpc > nine_index) && tmpc != space_index &&
                                 tmpc != positive_index) || nsigns > 1 ||
                     (number_begins && tmpc == space_index) ||
                     (tmpc == positive_index && j > 0 && prev_tmpc != space_index));
	  prev_tmpc = tmpc;
        }
      }
      break;
    case NumberFormat::CHAR4:
      break;
    }
    if (problem) {
      const std::string callerline = (caller.size() > 0) ? "  Called by: " + caller + "." : "";
      const std::string taskline = (task.size() > 0) ? "  Task: " + task + "." : "";
      rtErr("Malformed number on line " + std::to_string(iline + 1) + " of file " +
            tf.getFileName() + "." + callerline + taskline, "readNumberSeries");
    }
    iline++;
    llim += n_this_line;
  }

  // Real numbers and integers can all be read as number series
  iline = start_line;
  llim = 0;
  const int dot_position = width - decimal - 1 - (4 * (format == NumberFormat::SCIENTIFIC));
  const int e_position = width - 4;
  const double decimal_scale = pow(10.0, -decimal);
  while (llim < n_values) {
    const int hlim = std::min(llim + values_per_line, n_values);
    const int n_this_line = hlim - llim;
    const int readbase = tfr.line_limits[iline];
    int multiplier = 1;
    if (format_understood && fixed_resources_exceeded == false) {
      switch (format) {
      case NumberFormat::SCIENTIFIC:
        for (int i = 0; i < n_this_line; i++) {
          line_integral[i] = 0;
          line_decimal[i] = 0;
          line_exponent[i] = 0;
          line_negative[i] = 0;
        }

        // First do the exponential part
	for (int j = 0; j < n_this_line; j++) {
	  const int rjpw = readbase + ((j + 1) * width);
          line_exponent[j] = ((tfr.text[rjpw - 2] - zero_index) * 10) +
                             tfr.text[rjpw - 1] - zero_index;
          line_exponent[j] *= 1 - (2 * (tfr.text[rjpw - 3] == negative_index));
	}

	// Now do the decimal part
        for (int i = e_position - 1; i > dot_position; i--) {
          for (int j = 0; j < n_this_line; j++) {
            const int tmpc = tfr.text[readbase + (j * width) + i];
            line_decimal[j] += (tmpc - zero_index) * multiplier;
          }
          multiplier *= 10;
        }

	// Now do the integral part
	multiplier = 1;
        for (int i = dot_position - 1; i >= 0; i--) {
          for (int j = 0; j < n_this_line; j++) {
            const int tmpc = tfr.text[readbase + (j * width) + i];
            line_integral[j] += (tmpc >= zero_index && tmpc <= nine_index) *
                                (tmpc - zero_index) * multiplier;
            line_negative[j] += (tmpc == negative_index);
          }
          multiplier *= 10;
        }

	// Evaluate the components
        for (int i = llim; i < hlim; i++) {
          const int neg_flip = 1 - (2 * line_negative[i - llim]);
          line_integral[i - llim] *= neg_flip;
          line_decimal[i - llim] *= neg_flip;
          result[i].d = ((static_cast<double>(line_integral[i - llim]) +
                          (static_cast<double>(line_decimal[i - llim]) * decimal_scale)) *
                         pow(10.0, line_exponent[i - llim]));
        }

        break;
      case NumberFormat::STANDARD_REAL:
        for (int i = 0; i < n_this_line; i++) {
          line_integral[i] = 0;
          line_decimal[i] = 0;
          line_negative[i] = 0;
        }

	// First do the decimal part
        for (int i = width - 1; i > dot_position; i--) {
          for (int j = 0; j < n_this_line; j++) {
            const int tmpc = tfr.text[readbase + (j * width) + i];
            line_decimal[j] += (tmpc - zero_index) * multiplier;
          }
          multiplier *= 10;
        }

	// Now do the integral part
	multiplier = 1;
        for (int i = dot_position - 1; i >= 0; i--) {
          for (int j = 0; j < n_this_line; j++) {
            const int tmpc = tfr.text[readbase + (j * width) + i];
            line_integral[j] += (tmpc >= zero_index && tmpc <= nine_index) *
                                (tmpc - zero_index) * multiplier;
            line_negative[j] += (tmpc == negative_index);
          }
          multiplier *= 10;
        }

	// Add the two components
        for (int i = llim; i < hlim; i++) {
	  const int neg_flip = 1 - (2 * line_negative[i - llim]);
	  line_integral[i - llim] *= neg_flip;
	  line_decimal[i - llim] *= neg_flip;
          result[i].d = (static_cast<double>(line_integral[i - llim]) +
                         (static_cast<double>(line_decimal[i - llim]) * decimal_scale));
        }
        break;
      case NumberFormat::INTEGER:
        for (int i = 0; i < n_this_line; i++) {
          line_integral[i] = 0;
        }
        for (int i = width - 1; i >= 0; i--) {
          for (int j = 0; j < n_this_line; j++) {
            const int tmpc = tfr.text[readbase + (j * width) + i];
            line_integral[j] += (tmpc >= zero_index && tmpc <= nine_index) *
                                (tmpc - zero_index) * multiplier;
            line_integral[j] *= 1 - (2 * (tmpc == negative_index));
          }
          multiplier *= 10;
        }
        for (int i = llim; i < hlim; i++) {
          result[i].i = line_integral[i - llim];
        }
	break;
      case NumberFormat::LONG_LONG_INTEGER:
	{
          llint llmultiplier = 1LL;
          for (int i = 0; i < n_this_line; i++) {
            line_llintegral[i] = 0LL;
	    line_negative[i] = 0;
          }
          for (int i = width - 1; i >= 0; i--) {
            for (int j = 0; j < n_this_line; j++) {
              const int tmpc = tfr.text[readbase + (j * width) + i];
              line_llintegral[j] += static_cast<llint>((tmpc >= zero_index && tmpc <= nine_index) *
                                                       (tmpc - zero_index)) * llmultiplier;
              line_negative[j] += (tmpc == negative_index);
            }
            llmultiplier *= 10LL;
          }
          for (int i = llim; i < hlim; i++) {
            result[i].lli = line_llintegral[i - llim] *
                            static_cast<llint>(1 - (2 * line_negative[i - llim]));
          }
	}
	break;
      case NumberFormat::UNSIGNED_INTEGER:
	{
          ulint ulmultiplier = 1ul;
          for (int i = 0; i < n_this_line; i++) {
            line_uintegral[i] = 0ul;
          }
          for (int i = width - 1; i >= 0; i--) {
            for (int j = 0; j < n_this_line; j++) {
              const int tmpc = tfr.text[readbase + (j * width) + i];
              line_uintegral[j] += static_cast<ulint>((tmpc >= zero_index && tmpc <= nine_index) *
                                                       (tmpc - zero_index)) * ulmultiplier;
            }
            ulmultiplier *= 10ul;
          }
          for (int i = llim; i < hlim; i++) {
            result[i].ui = line_uintegral[i - llim];
          }
	}
	break;
      case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
	{
          ullint ullmultiplier = 1ull;
          for (int i = 0; i < n_this_line; i++) {
            line_ullintegral[i] = 0ull;
          }
          for (int i = width - 1; i >= 0; i--) {
            for (int j = 0; j < n_this_line; j++) {
              const int tmpc = tfr.text[readbase + (j * width) + i];
              line_ullintegral[j] += static_cast<ullint>((tmpc >= zero_index &&
                                                          tmpc <= nine_index) *
                                                         (tmpc - zero_index)) * ullmultiplier;
            }
            ullmultiplier *= 10ull;
          }
          for (int i = llim; i < hlim; i++) {
            result[i].ulli = line_ullintegral[i - llim];
          }
	}
	break;
      case NumberFormat::CHAR4:
	for (int i = 0; i < n_this_line; i++) {
          char4 tmpc4;
          tmpc4.x = tfr.text[readbase + 4*i];
          tmpc4.y = tfr.text[readbase + 4*i + 1];
          tmpc4.z = tfr.text[readbase + 4*i + 2];
          tmpc4.w = tfr.text[readbase + 4*i + 3];
	  result[llim + i].c4 = tmpc4;
	}
        break;
      }
    }
    else {
      for (int i = 0; i < n_this_line; i++) {
        for (int j = 0; j < width; j++) {
          buffer[j] = tfr.text[readbase + (i * width) + j];
        }
        switch (format) {
        case NumberFormat::SCIENTIFIC:
        case NumberFormat::STANDARD_REAL:
          result[llim + i].d = stod(buffer);
	  break;
        case NumberFormat::INTEGER:
          result[llim + i].i = stol(buffer);
	  break;
        case NumberFormat::LONG_LONG_INTEGER:
          result[llim + i].lli = stoll(buffer);
	  break;
        case NumberFormat::UNSIGNED_INTEGER:
          result[llim + i].ui = stoul(buffer);
	  break;
        case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
          result[llim + i].ulli = stoull(buffer);
	  break;
        case NumberFormat::CHAR4:
	  break;
	}
      }
    }
    iline++;
    llim += n_this_line;
  }

  return result;
}

//-------------------------------------------------------------------------------------------------
llint readIntegerValue(const char* number_text, const int start_index, const int number_length) {
  std::string buffer;
  buffer.resize(number_length);
  const int ilim = start_index + number_length;
  for (int i = start_index; i < ilim; i++) {
    buffer[i - start_index] = number_text[i];
  }
  return strtoll(buffer.c_str(), nullptr, 10);
}

//-------------------------------------------------------------------------------------------------
llint readIntegerValue(const std::string &number_text, const int start_index,
                       const int number_length) {
  std::string buffer = number_text.substr(start_index, number_length);
  return strtoll(buffer.c_str(), nullptr, 10);
}

//-------------------------------------------------------------------------------------------------
llint readIntegerValue(const TextFile &tf, const int line_index, const int start_index,
                       const int number_length) {
  std::string buffer;
  buffer.resize(number_length);

  // Check that the request makes sense
  if (line_index >= tf.getLineCount()) {
    rtErr("Request for integer data on line " + std::to_string(line_index) + " is invalid for a "
          "text file containining " + std::to_string(tf.getLineCount()) + " lines.",
          "readIntegerValue");
  }
  const int lstart = tf.getLineLimits(line_index);
  const int lend   = tf.getLineLimits(line_index + 1);
  if (lstart + start_index + number_length > lend) {
    rtErr("Request for integer data on line " + std::to_string(line_index) + " is invalid for a "
          "line containining " + std::to_string(lend - lstart) + " characters.",
          "readIntegerValue");
  }
  const char* number_ptr = tf.getTextPointer(tf.getLineLimits(line_index));
  if (verifyNumberFormat(number_ptr, NumberFormat::INTEGER, start_index, number_length) == false) {
    rtErr("Request for integer data on line " + std::to_string(line_index) + " encounters an "
          "invalid format.", "readIntegerValue");
  }
  const int ilim = start_index + number_length;
  for (int i = start_index; i < ilim; i++) {
    buffer[i - start_index] = number_ptr[i];
  }
  return strtoll(buffer.c_str(), nullptr, 10);
}

//-------------------------------------------------------------------------------------------------
double readRealValue(const char* number_text, const int start_index, const int number_length) {
  std::string buffer;
  buffer.resize(number_length);
  const int ilim = start_index + number_length;
  for (int i = start_index; i < ilim; i++) {
    buffer[i - start_index] = number_text[i];
  }
  return strtod(buffer.c_str(), nullptr);
}

//-------------------------------------------------------------------------------------------------
double readRealValue(const std::string &number_text, const int start_index,
                     const int number_length) {
  std::string buffer = number_text.substr(start_index, number_length);
  return strtod(buffer.c_str(), nullptr);
}

//-------------------------------------------------------------------------------------------------
double readRealValue(const TextFile &tf, const int line_index, const int start_index,
                     const int number_length) {
  std::string buffer;
  buffer.resize(number_length);

  // Check that the request makes sense
  if (line_index >= tf.getLineCount()) {
    rtErr("Request for integer data on line " + std::to_string(line_index) + " is invalid for a "
          "text file containining " + std::to_string(tf.getLineCount()) + " lines.",
          "readRealValue");
  }
  const int lstart = tf.getLineLimits(line_index);
  const int lend   = tf.getLineLimits(line_index + 1);
  if (lstart + start_index + number_length > lend) {
    rtErr("Request for integer data on line " + std::to_string(line_index) + " is invalid for a "
          "line containining " + std::to_string(lend - lstart) + " characters.",
          "readRealValue");
  }
  const char* number_ptr = tf.getTextPointer(tf.getLineLimits(line_index));
  if (verifyNumberFormat(number_ptr, NumberFormat::STANDARD_REAL, start_index,
                         number_length) == false &&
      verifyNumberFormat(number_ptr, NumberFormat::STANDARD_REAL, start_index,
                         number_length) == false) {
    rtErr("Request for real-valued data on line " + std::to_string(line_index) + " encounters an "
          "invalid format.", "readIntegerValue");
  }
  const int ilim = start_index + number_length;
  for (int i = start_index; i < ilim; i++) {
    buffer[i - start_index] = number_ptr[i];
  }
  return strtod(buffer.c_str(), nullptr);
}

//-------------------------------------------------------------------------------------------------
  
} // namespace parse
} // namespace stormm
