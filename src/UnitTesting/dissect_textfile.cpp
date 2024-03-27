#include "copyright.h"
#include "Math/statistical_enumerators.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "dissect_textfile.h"

namespace stormm {
namespace testing {

using parse::countDelimiters;
using parse::countWords;
using parse::digestText;
using parse::markGuardedText;
using parse::separateText;
using parse::TextFileReader;
using stmath::maxValue;
using stmath::mean;
using stmath::variance;
using stmath::VarianceMethod;

//-------------------------------------------------------------------------------------------------
double2 meanWordsPerLine(const TextFile &tf) {
  const TextFileReader tfr = tf.data();
  std::vector<int> wcounts(tfr.line_count);
  for (int i = 0; i < tfr.line_count; i++) {
    wcounts[i] = countWords(tfr.text, tfr.line_limits[i], tfr.line_lengths[i]);
  }
  const double2 result = { mean(wcounts), variance(wcounts, VarianceMethod::STANDARD_DEVIATION) };
  return result;
}
  
//-------------------------------------------------------------------------------------------------
int countInstances(const TextFile &tf, const std::string &sequence,
                   const std::vector<char> &separators, const std::vector<char> &delimiters,
                   const std::vector<TextGuard> &quote_marks,
                   const std::vector<TextGuard> &comment_marks) {
  const TextFileReader tfr = tf.data();

  // Begin by creating comment and quotation masks.
  const std::vector<bool> commented = markGuardedText(tfr, comment_marks, quote_marks);
  const std::vector<bool> quoted = markGuardedText(tfr, quote_marks, comment_marks);

  // Step through the entire text file, identifying separators not masked by quotations or comments
  // and marking them off for additional spaces.  This is similar to what happens with one stage of
  // &namelist parsing, but the code is not so significant that it warrants its own encapsulation.
  const int longest_line = maxValue(tfr.line_lengths, tfr.line_count);
  std::string buffer;
  buffer.reserve(longest_line);
  const size_t seqlen = sequence.size();
  const size_t nsep = separators.size();
  const size_t ndelim = delimiters.size();
  int instance_count = 0;
  std::vector<bool> buffer_quoted, buffer_commented;
  for (int i = 0; i < tfr.line_count; i++) {
    buffer = digestText(tfr.text, tfr.line_limits[i], tfr.line_lengths[i], separators.data(),
                        separators.size(), delimiters.data(), delimiters.size(), quoted, commented,
                        &buffer_quoted, &buffer_commented);

    // Search for the word of interest within the buffer.
    const size_t search_limit = buffer.size() - seqlen;
    bool word_found = false;
    int search_start = 0;
    while (search_start <= search_limit) {
      size_t k = 0;
      bool search_valid;
      do {
        const size_t ssk = search_start + k;
        search_valid = (k < seqlen && buffer[ssk] == sequence[k] && buffer_quoted[ssk] == false &&
                        buffer_commented[ssk] == false);
        if (search_valid) {
          k++;
        }
      } while (search_valid);
      if (k == seqlen) {
        
        // The word is found in this line.  However, it must be a unique word separated by spaces.
        if ((search_start == 0 || buffer[search_start - 1] == ' ') &&
            (search_start == search_limit || buffer[search_start + k] == ' ')) {
          instance_count++;
        }
      }
      search_start++;
    }
  }
  return instance_count;
}

//-------------------------------------------------------------------------------------------------
int countInstances(const TextFile &tf, const std::vector<std::string> &sequence,
                   const std::vector<char> &separators, const std::vector<char> &delimiters,
                   const std::vector<TextGuard> &quote_marks,
                   const std::vector<TextGuard> &comment_marks) {

  // Error if no search sequence is provided
  if (sequence.size() == 0) {
    rtErr("No search sequence was provided for searching text file " + tf.getFileName() + ".",
          "countInstances");
  }
  size_t total_seq_length = 0;
  const size_t seq_word_count = sequence.size();
  for (size_t i = 0; i < seq_word_count; i++) {
    total_seq_length += sequence[i].size();
  }
  total_seq_length += seq_word_count - 1;
  if (total_seq_length == 0) {
    rtErr("No searchable characters were provided for searching text file " + tf.getFileName() +
          ".", "countInstances");
  }
  const TextFileReader tfr = tf.data();

  // Begin by converting the TextFile to a string, then making comment and quotation masks.
  const std::string buffer = tf.toString(TextEnds::NEWLINE);
  const std::vector<bool> commented = markGuardedText(buffer, comment_marks, quote_marks);
  const std::vector<bool> quoted = markGuardedText(buffer, quote_marks, comment_marks);
  
  // Create a copy of the text file, placing white space for every line termination but otherwise
  // doing away with implicit carriage returns.
  std::vector<bool> sp_quoted, sp_commented;
  const std::string sp_buffer = digestText(buffer.data(), 0, buffer.size(), separators.data(),
                                           separators.size(), delimiters.data(), delimiters.size(),
                                           quoted, commented, &sp_quoted, &sp_commented);
  
  // The search must extend over all words in the sequence, in order
  const size_t spbuff_length = sp_buffer.size();
  const size_t search_limit = spbuff_length - total_seq_length;
  size_t search_start = 0;
  int instance_count = 0;
  while (search_start <= search_limit) {
    size_t i = 0;
    size_t k = 0;
    
    // The search can only be valid if the window is preceded by white space.
    bool search_valid = ((search_start == 0 || sp_buffer[search_start - 1] == ' ') &&
                         sp_buffer[search_start] != ' ' && sp_quoted[search_start] == false &&
                         sp_commented[search_start] == false);
    while (search_valid && i < seq_word_count) {
      
      // Step forward to the start of the next word
      while (search_start +  k < spbuff_length && sp_buffer[search_start + k] == ' ') {
        k++;
      }
      const size_t seqi_length = sequence[i].size();
      if (search_start + k <= spbuff_length - seqi_length) {
        size_t j = 0;
        size_t sp_jk = search_start + j + k;
        while (j < seqi_length && sp_buffer[sp_jk] == sequence[i][j] &&
               sp_quoted[sp_jk] == false && sp_commented[sp_jk] == false) {
          j++;
          sp_jk++;
        }
        if (j == seqi_length) {

          // The word was found.  Check that the following character is either the end of the
          // spaced buffer, commented, quoted, or a white space character.
          if (sp_jk == spbuff_length || sp_buffer[sp_jk] == ' ' || sp_commented[sp_jk] ||
              sp_quoted[sp_jk]) {
            i++;
            k += seqi_length;
          }
          else {
            search_valid = false;
          }
        }
        else {
          search_valid = false;
        }
      }
      else {
        search_valid = false;
      }
    }
    if (i == seq_word_count) {
      instance_count++;
    }
    search_start++;
  }
  return instance_count;
}

//-------------------------------------------------------------------------------------------------
TextFile grep(const TextFile &tf, const std::regex &criterion, const int before,
              const int after, const CaseSensitivity cased,
              const std::vector<TextGuard> &quote_marks,
              const std::vector<TextGuard> &comment_marks) {
  const TextFileReader tfr = tf.data();
  std::vector<bool> hits(tfr.line_count, false);
  const std::vector<bool> commented = markGuardedText(tfr, comment_marks, quote_marks);
  const std::vector<bool> quoted = markGuardedText(tfr, quote_marks, comment_marks);
  bool comments_active = false;
  bool quotes_active = false;
  for (size_t i = 0; i < tfr.line_limits[tfr.line_count]; i++) {
    comments_active = (comments_active || commented[i]);
    quotes_active = (quotes_active || quoted[i]);
  }
  for (int i = 0; i < tfr.line_count; i++) {
    std::string tf_line = tf.getLineAsString(i);
    if (comments_active || quotes_active) {
      std::string buffer;
      const size_t ln_length = tfr.line_lengths[i];
      const size_t ln_offset = tfr.line_limits[i];
      int n_char = 0;
      for (size_t j = 0; j < ln_length; j++) {
        if (commented[j + ln_offset] == false && quoted[j + ln_offset] == false) {
          n_char++;
        }
      }
      buffer.reserve(n_char);
      for (size_t j = 0; j < ln_length; j++) {
        if (commented[j + ln_offset] == false && quoted[j + ln_offset] == false) {
          buffer.push_back(tfr.text[j + ln_offset]);
        }
      }
      tf_line = buffer;
    }
    if (std::regex_match(tf_line, criterion)) {
      hits[i] = true;
      for (int j = std::max(i - before, 0); j < i; j++) {
        hits[j] = true;
      }
      for (int j = i + 1; j < std::min(i + 1 + after, tfr.line_count); j++) {
        hits[j] = true;
      }
    }
  }
  size_t result_length = 0;
  for (int i = 0; i < tfr.line_count; i++) {
    if (hits[i]) {
      result_length += static_cast<size_t>(tfr.line_lengths[i] + 1);
    }
  }
  std::string result(result_length, ' ');
  size_t rcon = 0;
  for (int i = 0; i < tfr.line_count; i++) {
    if (hits[i]) {
      const size_t llim = tfr.line_limits[i];
      const size_t hlim = tfr.line_limits[i + 1];
      for (size_t j = llim; j < hlim; j++) {
        result[rcon] = tfr.text[j];
        rcon++;
      }
      result[rcon] = '\n';
      rcon++;
    }
  }
  TextFile findings(result, TextOrigin::RAM);
  return findings;
}

//-------------------------------------------------------------------------------------------------
TextFile grep(const TextFile &tf, const std::string &criterion, const int before, const int after,
              const CaseSensitivity cased, const std::vector<TextGuard> &quote_marks,
              const std::vector<TextGuard> &comment_marks) {
  return grep(tf, std::regex(criterion), before, after, cased, quote_marks, comment_marks);
}

//-------------------------------------------------------------------------------------------------
TextFile grep(const TextFile &tf, const char* criterion, const int length, const int before,
              const int after, const CaseSensitivity cased,
              const std::vector<TextGuard> &quote_marks,
              const std::vector<TextGuard> &comment_marks) {
  const std::string str_criterion(criterion, length);
  return grep(tf, std::regex(str_criterion), before, after, cased, quote_marks, comment_marks);
}

} // namespace testing
} // namespace stormm
