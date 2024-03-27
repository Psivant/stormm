#include <limits.h>
#include <string>
#include <sys/ioctl.h>
#include <unistd.h>
#include "copyright.h"
#include "Reporting/display.h"
#include "Reporting/error_format.h"
#include "Reporting/summary_file.h"
#include "input_transcript.h"
#include "namelist_enumerators.h"
#include "nml_conformer.h"
#include "nml_dynamics.h"
#include "nml_ffmorph.h"
#include "nml_files.h"
#include "nml_minimize.h"
#include "nml_precision.h"
#include "nml_random.h"
#include "nml_report.h"
#include "nml_restraint.h"
#include "nml_solvent.h"

namespace stormm {
namespace namelist {

using display::horizontalRule;
using review::indentText;

//-------------------------------------------------------------------------------------------------
std::string commandLineAsString(const UserSettings &ui, int width) {
  std::string result;
  if (width > 0) {
    result += horizontalRule("//", "", width);
  }
  const std::vector<std::string> args = ui.getCommandLineArguments();
  const size_t nargs = args.size();
  std::string arg_cat;
  for (size_t i = 0; i < nargs; i++) {
    arg_cat += args[i];
    if (i < nargs - 1) {
      arg_cat += " ";
    }
  }
  result += indentText(arg_cat, 2, width);
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string prepareInputTranscript(const UserSettings &ui, const int width,
                                   const int max_repetitions,
                                   const std::vector<NamelistEmulator> &extra_nml) {
  FilesControls ficon = ui.getFilesNamelistInfo();
  const int actual_width = (width >= 0) ? width : ui.getReportNamelistInfo().getReportFileWidth();
  std::string result;
  
  // Write the &files namelist, already extracted, first.
  const NamelistIntroduction cheader = NamelistIntroduction::COMPACT_HEADER;
  result += ficon.getTranscript().printContents(actual_width, max_repetitions, cheader);

  // Write namelists associated with with processes that integrate capabilities
  if (ui.getConformerPresence()) {
    const ConformerControls& confcon = ui.getConformerNamelistInfo();
    result += confcon.getTranscript().printContents(actual_width, max_repetitions, cheader);
  }
  if (ui.getDynamicsPresence()) {
    const DynamicsControls& dynacon = ui.getDynamicsNamelistInfo();
    result += dynacon.getTranscript().printContents(actual_width, max_repetitions, cheader);
  }
  if (ui.getFFMorphPresence()) {
    const FFMorphControls& morphcon = ui.getFFMorphNamelistInfo();
    result += morphcon.getTranscript().printContents(actual_width, max_repetitions, cheader);
  }
  
  // Write additional namelists
  if (ui.getReportPresence()) {
    const ReportControls& repcon = ui.getReportNamelistInfo();
    result += repcon.getTranscript().printContents(actual_width, max_repetitions, cheader);
  }
  if (ui.getMinimizePresence()) {
    const MinimizeControls& mincon = ui.getMinimizeNamelistInfo();
    result += mincon.getTranscript().printContents(actual_width, max_repetitions, cheader);
  }
  if (ui.getSolventPresence()) {
    const SolventControls& gbcon = ui.getSolventNamelistInfo();
    result += gbcon.getTranscript().printContents(actual_width, max_repetitions, cheader);
  }
  if (ui.getRandomPresence()) {
    const RandomControls& rngcon = ui.getRandomNamelistInfo();
    result += rngcon.getTranscript().printContents(actual_width, max_repetitions, cheader);
  }
  if (ui.getPrecisionPresence()) {
    const PrecisionControls& preccon = ui.getPrecisionNamelistInfo();
    result += preccon.getTranscript().printContents(actual_width, max_repetitions, cheader);
  }

  // Write the series of all NMR restraint namelists
  const std::vector<RestraintControls>& rstcon = ui.getRestraintNamelistInfo();
  for (size_t i = 0; i < rstcon.size(); i++) {
    if (i == 0) {
      result += rstcon[i].getTranscript().printContents(actual_width, max_repetitions, cheader);
    }
    else {
      result += rstcon[i].getTranscript().printContents(actual_width, max_repetitions,
                                                        NamelistIntroduction::BLANK_LINE);
    }
  }

  // Write the content of additional namelists and return the result.
  for (size_t i = 0; i < extra_nml.size(); i++) {
    result += extra_nml[i].printContents(actual_width, max_repetitions, cheader);
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
void writeInputTranscript(const UserSettings &ui, const int width,
                          const int max_repetitions,
                          const std::vector<NamelistEmulator> &extra_nml) {
  FilesControls ficon = ui.getFilesNamelistInfo();
  const int actual_width = (width >= 0) ? width : ui.getReportNamelistInfo().getReportFileWidth();

  // If there is no transcript file specified, return immediately
  const std::string& xscript_name = ficon.getInputTranscriptFile();
  if (xscript_name.size() == 0) {
    return;
  }
  std::ofstream foutp = openOutputFile(xscript_name, ui.getPrintingPolicy(), "print a transcript "
                                       "of the user input");
  const std::string kw_text = prepareInputTranscript(ui, width, max_repetitions, extra_nml);
  foutp.write(kw_text.data(), kw_text.size());
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
void writeInputTranscript(const UserSettings &ui, const int max_repetitions,
                          const std::vector<NamelistEmulator> &extra_nml) {

  // Determine the file width based on the present terminal.
  struct winsize console_dims;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &console_dims);
  const int width = console_dims.ws_col - 1;
  writeInputTranscript(ui, width, max_repetitions, extra_nml);
}

} // namespace namelist
} // namespace stormm
