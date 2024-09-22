#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sys/time.h>
#include "copyright.h"
#include "Math/summation.h"
#include "Math/statistics.h"
#include "Parsing/parse.h"
#include "Parsing/tabulation.h"
#include "Reporting/error_format.h"
#include "stopwatch.h"

namespace stormm {
namespace testing {

using stmath::sum;
using parse::BorderFormat;
using parse::findStringInVector;
using parse::polyNumericVector;
using parse::NumberFormat;
using parse::printTable;
using parse::TableHeadingLine;

//-------------------------------------------------------------------------------------------------
StopWatch::StopWatch(const std::string &label_in) :
    category_count{1},
    initial_time{translateCurrentTime()},
    time_at_last_test{initial_time},
    category_start_times{std::vector<double>(16, initial_time)},
    category_total_times{std::vector<double>(16, 0.0)},
    category_squared_times{std::vector<double>(16, 0.0)},
    category_last_times{std::vector<double>(16, initial_time)},
    category_min_interval{std::vector<double>(16, 0.0)},
    category_max_interval{std::vector<double>(16, 0.0)},
    category_samples{std::vector<int>(16, 0)},
    category_names{std::vector<std::string>(16, "")},
    label{label_in}
{
  // The first section shall be called "Miscesllaneous." If no time intervals are assigned to
  // it, this auto-generated section will not even be reported.
  category_names[0] = "Miscellaneous";
  category_names.resize(category_count);
  
  // The capacities of various vectors are preset to be spacious.  Set the actual sizes back to 1.
  category_start_times.resize(category_count);
  category_total_times.resize(category_count);
  category_squared_times.resize(category_count);
  category_last_times.resize(category_count);
  category_min_interval.resize(category_count);
  category_max_interval.resize(category_count);
  category_samples.resize(category_count);
}

//-------------------------------------------------------------------------------------------------
double StopWatch::getTimeAtStart() const {
  return initial_time;
}

//-------------------------------------------------------------------------------------------------
double StopWatch::getTimeAtLastTest() const {
  return time_at_last_test;
}

//-------------------------------------------------------------------------------------------------
double StopWatch::getTimeSinceStart() const {
  return translateCurrentTime() - initial_time;
}

//-------------------------------------------------------------------------------------------------
double StopWatch::getTimeSinceLastTest() const {
  return translateCurrentTime() - time_at_last_test; 
}

//-------------------------------------------------------------------------------------------------
double StopWatch::getCategoryDuration(const int query_index) const {
  validateCategoryIndex(query_index, "getCategoryDuration");
  return category_total_times[query_index];
}

//-------------------------------------------------------------------------------------------------
double StopWatch::getCategoryDuration(const std::string &query_name) const {
  return getCategoryDuration(validateCategoryName(query_name, "getCategoryDuration"));
}

//-------------------------------------------------------------------------------------------------
double StopWatch::getCategoryAverageInterval(const int query_index) const {
  validateCategoryIndex(query_index, "getCategoryAverageInterval");
  return category_total_times[query_index] / static_cast<double>(category_samples[query_index]);
}

//-------------------------------------------------------------------------------------------------
double StopWatch::getCategoryAverageInterval(const std::string &query_name) const {
  return getCategoryAverageInterval(validateCategoryName(query_name,
                                                         "getCategoryAverageInterval"));
}

//-------------------------------------------------------------------------------------------------
double StopWatch::getCategoryMinimumTime(const int query_index) const {
  validateCategoryIndex(query_index, "getCategoryMinimumTime");
  return category_min_interval[query_index];
}
  
double StopWatch::getCategoryMinimumTime(const std::string &query_name) const {
  return getCategoryMinimumTime(validateCategoryName(query_name, "getCategoryMinimumTime"));
}

//-------------------------------------------------------------------------------------------------
double StopWatch::getCategoryMaximumTime(const int query_index) const {
  validateCategoryIndex(query_index, "getCategoryMaximumTime");
  return category_max_interval[query_index];
}

//-------------------------------------------------------------------------------------------------
double StopWatch::getCategoryMaximumTime(const std::string &query_name) const {
  return getCategoryMaximumTime(validateCategoryName(query_name, "getCategoryMaximumTime"));
}

//-------------------------------------------------------------------------------------------------
int StopWatch::getCategorySamples(const int query_index) const {
  validateCategoryIndex(query_index, "getCategorySamples");
  return category_samples[query_index];
}

//-------------------------------------------------------------------------------------------------
int StopWatch::getCategorySamples(const std::string &query_name) const {
  return getCategoryMaximumTime(validateCategoryName(query_name, "getCategorySamples"));
}

//-------------------------------------------------------------------------------------------------
std::string StopWatch::getCategoryName(const int query_index) const {
  validateCategoryIndex(query_index, "getCategoryName");
  return category_names[query_index];
}

//-------------------------------------------------------------------------------------------------
int StopWatch::getCategoryIndex(const std::string &query, const ExceptionResponse policy) const {
  const int query_index = findStringInVector(category_names, query);
  if (query_index == category_count) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("No category \"" + query + "\" was found.", "StopWatch", "getCategoryIndex");
    case ExceptionResponse::WARN:
      rtWarn("No category \"" + query + "\" was found.  An index of -1 will be returned and may "
             "cause problems later in the program.", "StopWatch", "getCategoryIndex");
      return -1;
    case ExceptionResponse::SILENT:
      return -1;
    }
  }
  return query_index;
}

//-------------------------------------------------------------------------------------------------
double StopWatch::getTotalDuration() const {
  return sum<double>(category_total_times);
}

//-------------------------------------------------------------------------------------------------
void StopWatch::assignTime(const int target, const double reference_time) {
  double current_time = translateCurrentTime();
  double delta = current_time - reference_time;
  category_total_times[target] += delta;
  category_squared_times[target] += delta * delta;
  category_last_times[target] = current_time;
  if (category_samples[target] == 0) {
    category_min_interval[target] = delta;
    category_max_interval[target] = delta;
  }
  else {
    category_max_interval[target] = std::max(category_max_interval[target], delta);
    category_min_interval[target] = std::min(category_min_interval[target], delta);
  }
  category_samples[target] += 1;
  time_at_last_test = current_time;
}

//-------------------------------------------------------------------------------------------------
void StopWatch::assignTime(const int target) {
  assignTime(target, time_at_last_test);
}

//-------------------------------------------------------------------------------------------------
void StopWatch::assignTime(const std::string &target) {
  assignTime(validateCategoryName(target, "assignTime"), time_at_last_test);
}

//-------------------------------------------------------------------------------------------------
void StopWatch::assignTime(const std::string &target, const double reference_time) {
  assignTime(validateCategoryName(target, "assignTime"), reference_time);
}

//-------------------------------------------------------------------------------------------------
int StopWatch::addCategory(const std::string &name) {
  const int cat_idx = findStringInVector(category_names, name);
  if (cat_idx < category_count) {
    return cat_idx;
  }
  category_count += 1;
  double current_time = translateCurrentTime();
  category_start_times.push_back(current_time);
  category_total_times.push_back(0.0);
  category_squared_times.push_back(0.0);
  category_last_times.push_back(current_time);
  category_min_interval.push_back(0.0);
  category_max_interval.push_back(0.0);
  category_samples.push_back(0);
  category_names.push_back(name);
  return category_count - 1;
}

//-------------------------------------------------------------------------------------------------
void StopWatch::printResults(const double precision) {

  // Craft the header depending on the overall scale of the execution time
  const double overall_time = sum<double>(category_total_times);
  int n_reported_categories = 0;
  for (int i = 0; i < category_count; i++) {
    n_reported_categories += (category_samples[i] > 0);
  }
  if (overall_time < 600.0) {
    printf("Timings data for %s (%d categories, total %.4lf seconds):\n\n", label.c_str(),
           n_reported_categories, overall_time);
  }
  else if (overall_time < 3600.0) {
    const int minutes = overall_time / 60.0;
    const double seconds = overall_time - (60.0 * static_cast<double>(minutes));
    printf("Timings data for %s (%d categories, total %d minutes, %.4lf seconds:\n\n",
           label.c_str(), n_reported_categories, minutes, seconds);
  }
  else {
    const int hours = overall_time / 3600.0;
    const int minutes = (overall_time - (3600.0 * static_cast<double>(hours))) / 60.0;
    const double seconds = overall_time - (3600.0 * static_cast<double>(hours)) -
                           (60.0 * static_cast<double>(minutes));
    printf("Timings data for %s (%d categories, total %d hours, %d minutes, %.4lf seconds):\n\n",
           label.c_str(), n_reported_categories, hours, minutes, seconds);
  }

  // Decide on the formatting for each column.  Also compute mean times and other statistics in
  // this loop over the data.
  int category_name_fmt = 8;
  int sample_count_fmt = 7;
  int total_width_fmt = 0;
  int mean_width_fmt = 0;
  int min_width_fmt = 0;
  int max_width_fmt = 0;
  int overall_dec_fmt = 4;

  // Erase unused categories
  std::vector<double> print_start_times;
  std::vector<double> print_total_times;
  std::vector<double> print_squared_times;
  std::vector<double> print_last_times;
  std::vector<double> print_min_interval;
  std::vector<double> print_max_interval;
  std::vector<int> print_samples;
  std::vector<std::string> print_names;
  for (int i = 0; i < category_count; i++) {
    if (category_samples[i] > 0) {
      print_start_times.push_back(category_start_times[i]);
      print_total_times.push_back(category_total_times[i]);
      print_squared_times.push_back(category_squared_times[i]);
      print_last_times.push_back(category_last_times[i]);
      print_min_interval.push_back(category_min_interval[i]);
      print_max_interval.push_back(category_max_interval[i]);
      print_samples.push_back(category_samples[i]);
      print_names.push_back(category_names[i]);
    }
  }
  const int print_count = print_names.size();
  std::vector<double> mean_ctg_times(print_count, 0.0);
  std::vector<double> std_ctg_times(print_count, 0.0);
  
  for (int i = 0; i < print_count; i++) {
    mean_ctg_times[i] = print_total_times[i] / static_cast<double>(print_samples[i]);
    std_ctg_times[i] = stmath::running_stdev(print_squared_times[i], print_total_times[i],
                                           print_samples[i]);
  }
  std::vector<NumberFormat> fmt_key(6, NumberFormat::STANDARD_REAL);
  fmt_key[0] = NumberFormat::INTEGER;
  printTable({"Category Name", "Samples", "Total Time, s", "Mean Time, s", "Standard Deviation",
              "Minimum Time, s", "Maximum Time, s"}, print_names,
             {polyNumericVector(print_samples), polyNumericVector(print_total_times),
              polyNumericVector(mean_ctg_times), polyNumericVector(std_ctg_times),
              polyNumericVector(print_min_interval), polyNumericVector(print_max_interval)},
             fmt_key, std::string(""), BorderFormat::FULL);
}

//-------------------------------------------------------------------------------------------------
double StopWatch::translateCurrentTime() const {
  timeval x;
  gettimeofday(&x, nullptr);
  return static_cast<double>(x.tv_sec) + (static_cast<double>(x.tv_usec) * 1.0e-6);
}

//-------------------------------------------------------------------------------------------------
void StopWatch::validateCategoryIndex(const int index,
                                      const std::string &referring_function) const {
  if (index >= category_count || index < 0) {
    rtErr("Invalid category index " +  std::to_string(index) + " (valid indices 0 to " +
	  std::to_string(category_count - 1) + ").", "StopWatch", referring_function.c_str());
  }
}

//-------------------------------------------------------------------------------------------------
int StopWatch::validateCategoryName(const std::string &query,
                                    const std::string &referring_function) const {
  const int query_index = findStringInVector(category_names, query);
  if (query_index == category_count) {
    rtErr("Invalid section name \"" + query + "\".", "StopWatch", referring_function.c_str());
  }
  return query_index;
}
  
} // namespace testing
} // namespace stormm
