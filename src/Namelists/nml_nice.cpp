#include "copyright.h"
#include "input.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "namelist_element.h"
#include "nml_nice.h"

namespace stormm {
namespace namelist {

using parse::minimalRealFormat;
using parse::NumberFormat;
using parse::TextOrigin;
using parse::uppercase;
using parse::verifyContents;

//-------------------------------------------------------------------------------------------------
NiceControls::NiceControls(const ExceptionResponse policy_in) :
    policy{policy_in},
    workday_effort_level{default_workday_gpu_effort},
    running_check_interval{default_nice_running_check_interval},
    resume_check_interval{default_nice_running_check_interval},
    stash_gpu_data{false},
    worktime_start{0},
    worktime_end{0},
    work_days{},
    working_day_mask{std::vector<bool>(7, false)},
    nml_transcript{"nice"}
{
  // Load in a blank namelist so that certain keywords will be present, as if this were the means
  // by which the data was loaded.
  std::string tfs("&nice\n&end\n");
  TextFile tf(tfs, TextOrigin::RAM);
  int start_line = 0;
  bool found;
  nml_transcript = niceInput(tf, &start_line, &found, ExceptionResponse::SILENT);
}

//-------------------------------------------------------------------------------------------------
NiceControls::NiceControls(const TextFile &tf, int *start_line, bool *found_nml,
                           const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    NiceControls(policy_in)
{
  NamelistEmulator t_nml = niceInput(tf, start_line, found_nml,policy_in, wrap);
  nml_transcript = t_nml;

  // Define the work week
  if (t_nml.getKeywordStatus("working_day") == InputStatus::MISSING) {
    const DayOfTheWeek week_start = translateDayOfTheWeek(t_nml.getStringValue("workweek_start"));
    const DayOfTheWeek week_stop = translateDayOfTheWeek(t_nml.getStringValue("workweek_end"));
    const int iweek_start = static_cast<int>(week_start);
    const int iweek_stop = static_cast<int>(week_stop);
    for (int i = iweek_start; i <= iweek_stop; i++) {
      const DayOfTheWeek iday = static_cast<DayOfTheWeek>(i);
      work_days.push_back(iday);
    }
  }
  else {
    const std::vector<std::string> tmp_work_days = t_nml.getAllStringValues("working_day");
    const size_t nw_days = tmp_work_days.size();
    work_days.reserve(nw_days);
    for (size_t i = 0; i < nw_days; i++) {
      work_days.push_back(translateDayOfTheWeek(tmp_work_days[i]));
    }
  }
  const int work_week_length = work_days.size();
  for (int i = 0; i < work_week_length; i++) {
    working_day_mask[static_cast<size_t>(work_days[i])] = true;
  }

  // Define the hours of official business
  setWorkTimeStart(t_nml.getStringValue("workday_start"));
  setWorkTimeEnd(t_nml.getStringValue("workday_end"));

  // Define how to manage the GPU jobs during work time
  setGpuEffortLevel(t_nml.getRealValue("max_workday_gpu"));
  setRunningCheckInterval(t_nml.getRealValue("running_watch"));
  setResumeCheckInterval(t_nml.getRealValue("resume_watch"));
  stash_gpu_data = t_nml.getBoolValue("stash_gpu_arrays");
}

//-------------------------------------------------------------------------------------------------
std::string NiceControls::getWorkdayStart() const {
  return translateTimeInteger(worktime_start);
}

//-------------------------------------------------------------------------------------------------
std::string NiceControls::getWorkdayEnd() const {
  return translateTimeInteger(worktime_end);
}

//-------------------------------------------------------------------------------------------------
std::string NiceControls::getWorkWeek() const {
  std::string result;
  const size_t ndays = work_days.size();
  for (size_t i = 0; i < ndays; i++) {
    result += getEnumerationName(work_days[i]);
    if (i < ndays - 1) {
      result += ", ";
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
bool NiceControls::isWorkTimeNow() const {
  time_t raw_time = time(nullptr);
  tm* current_time = localtime(&raw_time);
  const int minute_now = (current_time->tm_hour * 60) + current_time->tm_min;
  return (minute_now >= worktime_start && minute_now < worktime_end &&
          working_day_mask[current_time->tm_wday]);
}

//-------------------------------------------------------------------------------------------------
void NiceControls::setGpuEffortLevel(const double workday_effort_level_in) {
  if (validateGpuEffortLevel(workday_effort_level_in)) {;
    workday_effort_level = workday_effort_level_in;
  }
}

//-------------------------------------------------------------------------------------------------
void NiceControls::setRunningCheckInterval(const double running_check_interval_in) {
  if (validateRunningCheckInterval(running_check_interval_in)) {
    running_check_interval = running_check_interval_in;
  }
}

//-------------------------------------------------------------------------------------------------
void NiceControls::setResumeCheckInterval(const double resume_check_interval_in) {
  if (validateResumeCheckInterval(resume_check_interval_in)) {
    resume_check_interval = resume_check_interval_in;
  }
}

//-------------------------------------------------------------------------------------------------
void NiceControls::setWorkTimeStart(const std::string &work_time_start_in) {
  worktime_start = translateTimeString(work_time_start_in);
}

//-------------------------------------------------------------------------------------------------
void NiceControls::setWorkTimeEnd(const std::string &work_time_end_in) {
  worktime_end = translateTimeString(work_time_end_in);
}

//-------------------------------------------------------------------------------------------------
bool NiceControls::validateGpuEffortLevel(const double workday_effort_level_in) const {
  if (workday_effort_level_in >= 0.0 && workday_effort_level_in <= 1.0) {
    return true;
  }
  else {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A GPU effort level of " + minimalRealFormat(workday_effort_level_in, 1.0e-2, true) +
            " is invalid.", "NiceControls", "validateGpuEffortLevel");
    case ExceptionResponse::WARN:
      rtWarn("A GPU effort level of " + minimalRealFormat(workday_effort_level_in, 1.0e-2, true) +
             " is invalid.  The current value of " +
             minimalRealFormat(workday_effort_level_in, 1.0e-2, true) + " will be maintained.",
             "NiceControls", "validateGpuEffortLevel");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    return false;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool NiceControls::validateRunningCheckInterval(const double running_check_interval_in) const {
  if (running_check_interval_in >= minimum_running_check_interval &&
      running_check_interval_in <= maximum_running_check_interval) {
    return true;
  }
  else {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("An interval of " + minimalRealFormat(running_check_interval_in, 1.0e-1, true) +
            " is not permitted for checking for competing GPU usage.", "NiceControls",
            "validateRunningCheckInterval");
    case ExceptionResponse::WARN:
      rtWarn("An interval of " + minimalRealFormat(running_check_interval_in, 1.0e-1, true) +
            " is not permitted for checking for competing GPU usage.  The current value of " +
             minimalRealFormat(running_check_interval, 1.0e-1, true) + " will be maintained.",
             "NiceControls", "validateRunningCheckInterval");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    return false;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool NiceControls::validateResumeCheckInterval(const double resume_check_interval_in) const {
  if (resume_check_interval_in >= minimum_resume_check_interval &&
      resume_check_interval_in <= maximum_resume_check_interval) {
    return true;
  }
  else {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("An interval of " + minimalRealFormat(resume_check_interval_in, 1.0e-1, true) +
            " is not permitted for checking for competing GPU usage.", "NiceControls",
            "validateResumeCheckInterval");
    case ExceptionResponse::WARN:
      rtWarn("An interval of " + minimalRealFormat(resume_check_interval_in, 1.0e-1, true) +
            " is not permitted for checking for competing GPU usage.  The current value of " +
             minimalRealFormat(resume_check_interval, 1.0e-1, true) + " will be maintained.",
             "NiceControls", "validateResumeCheckInterval");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    return false;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int NiceControls::translateTimeString(const std::string &time_statement) const {
  std::string time_nocase = uppercase(time_statement);

  // Scan for an 'AM' or 'PM' at the end of the string
  size_t total_length = time_nocase.size();
  bool has_am = false;
  bool has_pm = false;
  bool problem = false;
  if (total_length >= 2 && time_nocase[total_length - 1] == 'M') {
    has_am = (time_nocase[total_length - 2] == 'A');
    has_pm = (time_nocase[total_length - 2] == 'P');
    if (total_length == 2) {
      problem = true;
    }
    total_length -= 2;
    time_nocase.resize(total_length);
  }

  // Scan for a ':' character
  std::string time_hours, time_minutes;
  bool has_colon = false;
  for (size_t i = 0; i < total_length; i++) {
    if (time_nocase[i] == ':') {
      has_colon = true;
      time_hours = time_nocase.substr(0, i);
      if (i < total_length - 1) {
        time_minutes = time_nocase.substr(i + 1, total_length);
      }
    }
  }
  if (has_colon == false) {
    time_hours = time_nocase;
  }
  problem = (problem || (time_hours.size() == 0));
  int hours, minutes;
  if (verifyContents(time_hours, NumberFormat::INTEGER)) {
    hours = stoi(time_hours);
  }
  else {
    problem = true;
  }
  if (time_minutes.size() > 0) {
    if (verifyContents(time_minutes, NumberFormat::INTEGER)) {
      minutes = stoi(time_minutes);
    }
    else {
      problem = false;
    }
  }
  else {
    minutes = 0;
  }
  if (problem) {
    rtErr("Invalid time \"" + time_statement + "\".", "NiceControls", "translateTimeString");
  }
  if (has_pm && hours < 12) {
    hours += 12;
  }
  return (hours * 60) + minutes;
}

//-------------------------------------------------------------------------------------------------
std::string NiceControls::translateTimeInteger(const int time_value) const {
  const int hours = time_value / 60;
  const int minutes = time_value - (hours * 60);
  std::string result;
  if (hours < 12) {
    if (minutes < 10) {
      result = std::to_string(hours) + ":0" + std::to_string(minutes) + " AM";
    }
    else {
      result = std::to_string(hours) + ":" + std::to_string(minutes) + " AM";
    }
  }
  else {
    if (hours > 12) {
      result = std::to_string(hours - 12);
    }
    else {
      result = std::to_string(hours);
    }
    if (minutes < 10) {
      result += ":0" + std::to_string(minutes) + " PM";
    }
    else {
      result += ":" + std::to_string(minutes) + " PM";
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator niceInput(const TextFile &tf, int *start_line, bool *found,
                           const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("nice", CaseSensitivity::AUTOMATIC, policy, "Collects job control and "
                         "priority settings to facilitate running STORMM processes on desktop "
                         "computing resources without interfering in the day-to-day work of the "
                         "user or a friend.");

  // Define the working day
  t_nml.addKeyword("workday_start", NamelistType::STRING,
                   std::to_string(default_workday_start_hour));
  t_nml.addKeyword("workday_end", NamelistType::STRING,
                   std::to_string(default_workday_end_hour));
  t_nml.addKeyword("workweek_start", NamelistType::STRING, std::string(default_workweek_start));
  t_nml.addKeyword("workweek_end", NamelistType::STRING, std::string(default_workweek_end));
  t_nml.addKeyword("working_day", NamelistType::STRING, std::string(""), DefaultIsObligatory::NO,
                   InputRepeats::YES);
  t_nml.addHelp("workday_start", "The beginning of the recognized workday.  If given as a single "
                "integer, the input will be converted to an hour in the range 0 to 23 (inclusive, "
                "this is \"military time\").  If a ':' character is detected, the string will be "
                "interpreted as hours and minutes.  A suffix of 'AM' or 'PM' (case-insensitive) "
                "will adjust an hour in the range 0 to 11 (inclusive) to morning or afternoon, "
                "respectively.");
  t_nml.addHelp("workday_end", "The end of the recognized workday.  If given as a single "
                "integer, the input will be converted to an hour in the range 0 to 23 (inclusive, "
                "this is \"military time\").  If a ':' character is detected, the string will be "
                "interpreted as hours and minutes.  A suffix of 'AM' or 'PM' (case-insensitive) "
                "will adjust an hour in the range 0 to 11 (inclusive) to morning or afternoon, "
                "respectively.");
  t_nml.addHelp("workweek_start", "The first day of the recognized working week.  Values of the "
                "day's full name, three letter abbreviation, or one to two letter abbreviation (S "
                "M T W Th F Sa) are all accepted, without case sensitivity.");
  t_nml.addHelp("workweek_end", "The last day of the recognized working week.  Values of the "
                "day's full name, three letter abbreviation, or one to two letter abbreviation (S "
                "M T W Th F Sa) are all accepted, without case sensitivity.");
  t_nml.addHelp("working_day", "Mark one day as a working day.  This can be useful when the "
                "\"work week\" does not come as a series of consecutive days.  Specifying one or "
                "more working days in this manner will cause the workweek_start and workweek_end "
                "keywords to be ignored.");

  // Define the GPU job management
  t_nml.addKeyword("max_workday_gpu", NamelistType::REAL,
                   std::to_string(default_workday_gpu_effort));
  t_nml.addHelp("max_workday_gpu", "Provide a value between 0 (no GPU usage will be permitted "
                "during working hours) and 1.0 (the GPU will be permitted to run at 100% during "
                "worktime hours).  The default is chosen to provide a balance between "
                "productivity and mitigation of fan noise.");
  t_nml.addKeyword("running_watch", NamelistType::REAL,
                   std::to_string(default_nice_running_check_interval));
  t_nml.addHelp("running_watch", "The interval, in seconds, at which to check for other activity "
                "while the STORMM process is running at whatever permitted usage of the GPU "
                "resources (100% outside of working hours, or at a user-specified effort level "
                "during the work day).");
  t_nml.addKeyword("resume_watch", NamelistType::REAL,
                   std::to_string(default_nice_resume_check_interval));
  t_nml.addHelp("resume_watch", "The interval, in seconds, at which to check for other activity "
                "if the STORMM process is ever paused to make room for a competing GPU job.");
  t_nml.addKeyword("stash_gpu_arrays", NamelistType::BOOLEAN);
  t_nml.addHelp("stash_gpu_arrays", "If specified, the GPU allocations of certain large objects "
                "will be stashed in temporary arrays on the CPU whenever a GPU job is temporarily "
                "paused.  This will help make room for other processes to take full advantage of "
                "the GPU when they exert priority.");
  
  // Find the namelist in the input, populate the emulator, and update the starting line for future
  // namelist reading.
  *start_line = readNamelist(tf, &t_nml, *start_line, wrap, tf.getLineCount(), found);
  return t_nml;
}

} // namespace namelist
} // namespace stormm
