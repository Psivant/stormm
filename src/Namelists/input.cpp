#include <string>
#include <vector>
#include "copyright.h"
#include "Parsing/parse.h"
#include "Parsing/textfile.h"
#include "Math/vector_ops.h"
#include "Reporting/error_format.h"
#include "input.h"

namespace stormm {
namespace namelist {

using parse::strncmpCased;
using parse::LineSpan;
using parse::TextFileReader;
using stmath::maxValue;

//-------------------------------------------------------------------------------------------------
std::string padNamelistTuples(const char* text, const std::vector<bool> &quoted,
                              std::vector<bool> &commented) {
  const int n_char = quoted.size();
  if (quoted.size() != commented.size()) {
    rtErr("Quotation and comment marking vectors must be the same size (" +
          std::to_string(quoted.size()) + " and " + std::to_string(commented.size()) +
          " were supplied).", "padNamelistTuples");
  }
  int n_visual_aids = 0;
  for (int i = 0; i < n_char; i++) {
    n_visual_aids += (quoted[i] == false && commented[i] == false &&
                      (text[i] == '=' || text[i] == ',' || text[i] == ';'));
  }
  std::string result(n_char + (2 * n_visual_aids), ' ');
  int j = 0;
  for (int i = 0; i < n_char; i++) {
    if (text[i] == '=' || text[i] == ',' || text[i] == ';') {
      result[j] = ' ';
      result[j + 1] = text[i];
      result[j + 2] = ' ';
      j += 3;
    }
    else {
      result[j] = text[i];
      j++;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> pullNamelist(const TextFile &tf, const NamelistEmulator &nml,
                                      const int start_line, const WrapTextSearch wrap,
                                      int *end_line) {
  
  // Make the scan, seeking out the title (plus leading ampersand) as a keyword
  const std::string title_key = "&" + nml.getTitle();
  const size_t n_title_char = title_key.size();
  if (n_title_char == 1) {
    rtErr("Cannot search for a namelist with no title in " + tf.getFileName() + ".",
          "pullNamelist");
  }
  int namelist_start = -1;
  const TextFileReader tfr = tf.data();
  const std::vector<TextGuard> comment_marks = { TextGuard("#"), TextGuard("!"), TextGuard("%"),
                                                 TextGuard("//"),
                                                 TextGuard("/*", "*/", LineSpan::MULTIPLE) };
  const std::vector<TextGuard> quote_marks = { TextGuard("\"", "\"", LineSpan::MULTIPLE),
                                               TextGuard("'", "'") };
  const std::vector<bool> commented = markGuardedText(tfr, comment_marks, quote_marks);
  const std::vector<bool> quoted = markGuardedText(tfr, quote_marks, comment_marks);
  const int max_line_length = maxValue(tfr.line_lengths, tfr.line_count);
  std::vector<char> buffer(max_line_length + 1, '\0');
  const int actual_end = (*end_line == -1 || *end_line > tfr.line_count) ?
                         tfr.line_count : *end_line;
  const CaseSensitivity case_sensitivity = nml.getCaseSensitivity();
  for (int i = start_line; i < actual_end; i++) {
    if (tfr.line_limits[i + 1] - tfr.line_limits[i] < n_title_char) {
      continue;
    }
    for (size_t j = tfr.line_limits[i]; j <= tfr.line_limits[i + 1] - n_title_char; j++) {

      // Check for the namelist title
      if (commented[j] == false && quoted[j] == false && tfr.text[j] == '&') {
        bool found = true;
        switch (case_sensitivity) {
        case CaseSensitivity::YES:
        case CaseSensitivity::NO:
          found = strncmpCased(title_key, &tfr.text[j], case_sensitivity);
          break;
        case CaseSensitivity::AUTOMATIC:
          found = strncmpCased(title_key, &tfr.text[j], CaseSensitivity::NO);
          break;
        }
        if (found) {

          // The namelist title has been found.  Process lines, first tweaking them to pad the
          // special tuple delimiters with whitespace (so that they will be identified as
          // individual words)
          std::vector<std::string> result;
          size_t start_char = j;
          int next_line = i + 1;
          bool end_found = false;
          std::string accumulated_nml;
          std::vector<bool> accumulated_comment_mask;
          std::vector<bool> accumulated_quote_mask;
          do {

            // Separate out delimiters for tuples and eliminate visual aids '=', ',' and ';'.
            int n_enclosures = 0;
            for (size_t k = start_char; k < tfr.line_limits[next_line]; k++) {
              n_enclosures += (quoted[k] == false && commented[k] == false &&
                               (tfr.text[k] == '[' || tfr.text[k] == '(' || tfr.text[k] == '{' ||
                                tfr.text[k] == ']' || tfr.text[k] == ')' || tfr.text[k] == '}'));
            }
            const int n_pad_char = tfr.line_limits[next_line] - start_char + (2 * n_enclosures) +
                                   1;
            std::vector<bool> line_quoted(n_pad_char, false);
            std::vector<bool> line_commented(n_pad_char, false);
            std::string tweaked_line(n_pad_char, ' ');
            int m = 0;
            for (int k = start_char; k < tfr.line_limits[next_line]; k++) {
              if (quoted[k] == false && commented[k] == false &&
                  (tfr.text[k] == '=' || tfr.text[k] == ',' || tfr.text[k] == ';')) {
                tweaked_line[m] = ' ';
                line_quoted[m] = false;
                line_commented[m] = false;
                m++;
              }
              else if (quoted[k] == false && commented[k] == false &&
                       (tfr.text[k] == '[' || tfr.text[k] == ']' || tfr.text[k] == '(' ||
                        tfr.text[k] == ')' || tfr.text[k] == '{' || tfr.text[k] == '}')) {
                for (int kp = 0; kp < 3; kp++) {
                  tweaked_line[m + kp] = ' ';
                  line_quoted[m + kp] = false;
                  line_commented[m + kp] = false;
                }
                tweaked_line[m + 1] = tfr.text[k];
                m += 3;
              }
              else {
                tweaked_line[m] = tfr.text[k];
                line_quoted[m] = quoted[k];
                line_commented[m] = commented[k];
                m++;
              }
            }

            // Extend the line accumulant by one character of white space.
            tweaked_line[m] = ' ';

            // Determine what to add to the quotation and comment masks.
            if (next_line < tfr.line_count) {
              line_quoted[m] = quoted[tfr.line_limits[next_line]];
              line_commented[m] = commented[tfr.line_limits[next_line]];
            }
            else {

              // If the file ends in an open quotation, it won't be part of a namelist, so this
              // editing will not interfere with normal operations of the comment and quotation
              // parsers.
              line_quoted[m] = false;
              line_commented[m] = false;
            }

            // Scan the refined line for a namelist terminating card.  This must be done manually
            // as multi-line comments could otherwise mess up separateText().  If a terminating
            // card is found, cut the line there, indicate that the end has been found, and
            // include only the part of the line up through the terminating card in the result.
            for (int k = 0; k < n_pad_char; k++) {
              if (line_commented[k] || line_quoted[k]) {
                continue;
              }
              if (tweaked_line[k] == '/' && (k == 0 || tweaked_line[k - 1] == ' ') &&
                  (k == n_pad_char - 1 || tweaked_line[k + 1] == ' ')) {
                tweaked_line = tweaked_line.substr(0, k + 1);
                end_found = true;
                break;
              }
              if (tweaked_line[k] == '&' && k < n_pad_char - 3) {
                bool card_clear = true;
                for (int m = 0; m < 4; m++) {
                  card_clear = (card_clear && line_commented[k + m] == false &&
                                line_quoted[k + m] == false);
                }
                if (card_clear) {
                  bool end_seems_found = false;
                  switch (case_sensitivity) {
                  case CaseSensitivity::YES:
                  case CaseSensitivity::NO:
                    if (strncmpCased(&tweaked_line[k], "&end", 4, case_sensitivity)) {
                      end_seems_found = true;
                    }
                  case CaseSensitivity::AUTOMATIC:
                    if (strncmpCased(&tweaked_line[k], "&end", 4, CaseSensitivity::NO)) {
                      end_seems_found = true;
                    }
                  }
                  if (end_seems_found) {

                    // Check that this &end card occurs as a separate word
                    if ((k > 0 && tweaked_line[k - 1] != ' ') ||
                        (k < n_pad_char - 4 && tweaked_line[k + 4] != ' ')) {
                      rtWarn("What appears to be an &end card is concatenated to additional text "
                             "on line " + std::to_string(next_line) + " of input from file " +
                             tf.getFileName() + ".  This may indicate a mistake in the input.",
                             "pullNamelist");
                    }
                    else {
                      tweaked_line = tweaked_line.substr(0, k + 4);
                      end_found = true;
                      break;
                    }
                  }
                }
              }
            }
            
            // Add this string, with visual aids removed and enclosures separated, to a
            // growing string containing the entire namelist.  The line itself gets one character
            // of added white space, and the comment / quotation masks are extended with the
            // status of the final character on the original line.
            accumulated_nml += tweaked_line;
            const int n_added_char = tweaked_line.size();
            for (int k = 0; k < n_added_char; k++) {
              accumulated_comment_mask.push_back(line_commented[k]);
              accumulated_quote_mask.push_back(line_quoted[k]);
            }
            start_char = tfr.line_limits[next_line];
            next_line++;
          } while (next_line <= tfr.line_count && end_found == false);

          // Parse words from the C++ string and masks produced by the above loop
          result = separateText(accumulated_nml, accumulated_comment_mask, accumulated_quote_mask,
                                quote_marks);

          // Check that the namelist was properly terminated
          if (end_found == false) {
            rtWarn("Unterminated namelist.  Input will be inferred from the content at hand.  "
                   "No guarantees as to the veracity of program execution are implied.",
                   "pullNamelist");
          }

          // Return this result: it was the first instance of the namelist
          *end_line = next_line;
          return result;
        }
      }
    }
  }

  // If wrapping is requested, there is one more opportunity to find a namelist.
  switch (wrap) {
  case WrapTextSearch::NO:
    {
      // Return an empty container to indicate that no namelist could be pulled
      *end_line = tf.getLineCount();
      std::vector<std::string> blank;
      return blank;
    }
    break;
  case WrapTextSearch::YES:
    {
      // The start line must be passed as a non-const variable, and the final line searched
      // must be passed back up.  This may leave the last line set to a value earlier than
      // where it started if wrapping is turned on.
      *end_line = start_line;
      return pullNamelist(tf, nml, 0, WrapTextSearch::NO, end_line);
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int readNamelist(const TextFile &tf, NamelistEmulator *nml, const int start_line,
                 WrapTextSearch wrap, const int end_line, bool *found)
{
  // Check that the namelist begins in its default state
  for (size_t i = 0; i < nml->getKeywordCount(); i++) {
    
    // Funny to have to get the keyword by index and then go back in to have the namelist search
    // for that keyword's index, but it's rare to have to run through the namelist's keywords
    // one by one, not knowing what they actually are.  Developers are not expected to access
    // namelist keywords by their indices in a list.
    const int n_entries = nml->getKeywordEntries(nml->getKeyword(i));
    if (n_entries > 1 && nml->getKeywordStatus(nml->getKeyword(i)) != InputStatus::DEFAULT) {
      rtErr("Namelist \"" + nml->getTitle() + "\" seems to already contain user-specified data.  "
            "There are " + std::to_string(n_entries) + " entries for keyword " +
            nml->getKeyword(i) + ".", "readNamelist");
    }
  }

  // Get the list of words for this namelist
  int nc_end_line = end_line;
  std::vector<std::string> nml_words = pullNamelist(tf, nml->getTitle(), start_line, wrap,
                                                    &nc_end_line);
  
  // Scan through the list, trying to identify keywords followed by values.  The first word of
  // any namelist's word vector will be the title card of the namelist itself, and the last word
  // will be "&end" or "/".  Neither needs interpretation.
  const int word_count = nml_words.size();
  if (found != nullptr) {
    *found = (word_count > 0);
  }
  const int nparam = nml->getKeywordCount();
  for (int i = 1; i < word_count - 2; i++) {

    // Check for the left side of an enclosure: if that's found, then the keyword needs to
    // correspond to a namelist STRUCT-type element.  Otherwise it must be an INTEGER, REAL,
    // or STRING.  The traps in code called by assignElement() will handle any errors.
    if (nml_words[i + 1] == "[" || nml_words[i + 1] == "(" || nml_words[i + 1] == "{") {
      if (nml->getKeywordKind(nml_words[i]) != NamelistType::STRUCT) {
        std::string message = "In namelist " + nml->getTitle() + ", keyword " + nml_words[i] +
                              " is followed by enclosure \"" + nml_words[i + 1] + "\" but is not "
                              "associated with STRUCT data.";
        switch (nml->getPolicy()) {
        case ExceptionResponse::DIE:
          rtErr(message, "readNamelist"); 
        case ExceptionResponse::WARN:
          rtWarn(message, "readNamelist"); 
        case ExceptionResponse::SILENT:
          break;
        }
      }
      int j = i + 2;
      while (nml_words[j] != "]" && nml_words[j] != ")" && nml_words[j] != "}" &&
             j < word_count - 1) {
        j += nml->assignElement(nml_words[i], nml_words[j], nml_words[j + 1]);
        j++;
      }
      if (nml_words[j] != "]" && nml_words[j] != ")" && nml_words[j] != "}") {
        std::string series;
        for (int k = i + 1; k < j; k++) {
          series += nml_words[k];
          if (k < j - 1) {
            series += " ";
          }
        }
        rtErr("No closure was found to the series for STRUCT namelist element " + nml_words[i] +
              ": " + series, "readNamelist");
      }
      nml->triggerResizeBuffer(nml_words[i]);
      i = j;
    }
    else {
      i += nml->assignElement(nml_words[i], nml_words[i + 1]);
    }
  }
  return nc_end_line;
}

} // namspace namelist
} // namespace stormm
