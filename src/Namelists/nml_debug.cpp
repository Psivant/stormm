#include "copyright.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "namelist_common.h"
#include "namelist_element.h"
#include "nml_debug.h"

namespace stormm {
namespace namelist {

using constants::CaseSensitivity;
using parse::NumberFormat;
using parse::realToString;
using parse::strcmpCased;
using parse::TextOrigin;

//-------------------------------------------------------------------------------------------------
DebugControls::DebugControls(const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    policy{policy_in}, check_forces{false}, report_large_forces{false},
    check_neighbor_list{false}, track_momentum_purge{false}, enforce_sanity{false},
    maximum_reports{default_max_anomaly_reports},
    inspection_interval{default_inspection_interval},
    large_force_threshold{default_large_force_threshold},
    high_speed_threshold{default_high_speed_threshold},
    bond_strain_sanity{default_bond_max_strain},
    angle_strain_sanity{default_angle_max_strain},
    varname{default_debug_variable_prefix},
    nml_transcript{"debug"}
{
  // Load in a blank namelist so that certain keywords will be present, as if this were the means
  // by which the data was loaded.
  std::string tfs("&debug\n&end\n");
  TextFile tf(tfs, TextOrigin::RAM);
  int start_line = 0;
  bool found;
  nml_transcript = debugInput(tf, &start_line, &found, ExceptionResponse::SILENT);
}

//-------------------------------------------------------------------------------------------------
DebugControls::DebugControls(const TextFile &tf, int *start_line, bool *found_nml,
                             const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    DebugControls(policy_in, wrap)
{
  NamelistEmulator t_nml = debugInput(tf, start_line, found_nml, policy, wrap);
  nml_transcript = t_nml;

  // There is no common keyword data for the debugging namelist.  Take in namelist-specific
  // keyword data.
  check_forces = t_nml.getBoolValue("forces");
  report_large_forces = t_nml.getBoolValue("force_trigger");
  check_neighbor_list = t_nml.getBoolValue("ngbr_placement");
  track_momentum_purge = t_nml.getBoolValue("track_purge");
  enforce_sanity = t_nml.getBoolValue("enforce_sanity");
  t_nml.assignVariable(&maximum_reports, "max_reports");
  t_nml.assignVariable(&inspection_interval, "interval_trigger");
  t_nml.assignVariable(&large_force_threshold, "force_threshold");
  t_nml.assignVariable(&high_speed_threshold, "speed_threshold");
  t_nml.assignVariable(&bond_strain_sanity, "bond_sane");
  t_nml.assignVariable(&angle_strain_sanity, "angle_sane");
}

//-------------------------------------------------------------------------------------------------
bool DebugControls::checkForces() const {
  return check_forces;
}

//-------------------------------------------------------------------------------------------------
bool DebugControls::reportLargeForces() const {
  return report_large_forces;
}

//-------------------------------------------------------------------------------------------------
bool DebugControls::checkNeighborListComp() const {
  return check_neighbor_list;
}

//-------------------------------------------------------------------------------------------------
bool DebugControls::trackMomentumPurge() const {
  return track_momentum_purge;
}

//-------------------------------------------------------------------------------------------------
int DebugControls::getMaximumReports() const {
  return maximum_reports;
}

//-------------------------------------------------------------------------------------------------
int DebugControls::getInspectionInterval() const {
  return inspection_interval;
}

//-------------------------------------------------------------------------------------------------
double DebugControls::getLargeForceThreshold() const {
  return large_force_threshold;
}

//-------------------------------------------------------------------------------------------------
double DebugControls::getHighSpeedThreshold() const {
  return high_speed_threshold;
}

//-------------------------------------------------------------------------------------------------
bool DebugControls::runSanityChecks() const {
  return enforce_sanity;
}

//-------------------------------------------------------------------------------------------------
double DebugControls::getBondStrainTolerance() const {
  return bond_strain_sanity;
}

//-------------------------------------------------------------------------------------------------
double DebugControls::getAngleStrainTolerance() const {
  return angle_strain_sanity;
}

//-------------------------------------------------------------------------------------------------
const std::string& DebugControls::getVariableBase() const {
  return varname;
}

//-------------------------------------------------------------------------------------------------
void DebugControls::setForceCheck(const bool setting_in) {
  check_forces = setting_in;
}

//-------------------------------------------------------------------------------------------------
void DebugControls::setLargeForceReport(const bool setting_in) {
  check_forces = setting_in;
}

//-------------------------------------------------------------------------------------------------
void DebugControls::setNeighborListCompCheck(const bool setting_in) {
  check_neighbor_list = setting_in;
}

//-------------------------------------------------------------------------------------------------
void DebugControls::setMomentumPurgeTracking(const bool setting_in) {
  track_momentum_purge = setting_in;
}

//-------------------------------------------------------------------------------------------------
void DebugControls::setMaximumReports(const int maximum_reports_in) {
  maximum_reports = maximum_reports_in;
}

//-------------------------------------------------------------------------------------------------
void DebugControls::setInspectionInterval(const int interval_in) {
  inspection_interval = interval_in;
}

//-------------------------------------------------------------------------------------------------
void DebugControls::setLargeForceThreshold(const double threshold_in) {
  large_force_threshold = threshold_in;
}

//-------------------------------------------------------------------------------------------------
void DebugControls::setHighSpeedThreshold(const double threshold_in) {
  high_speed_threshold = threshold_in;
}

//-------------------------------------------------------------------------------------------------
void DebugControls::setRunSanityChecks(const bool setting_in) {
  enforce_sanity = setting_in;
}

//-------------------------------------------------------------------------------------------------
void DebugControls::setVariableBase(const std::string &varname_in) {
  varname = varname_in;
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator debugInput(const TextFile &tf, int *start_line, bool *found,
                            const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("debug", CaseSensitivity::AUTOMATIC, policy, "Wraps directives needed "
                         "to carry out debugging of molecular simulations");

  // Add namelist-specific keywords
  t_nml.addKeyword("forces", NamelistType::BOOLEAN);
  t_nml.addKeyword("force_trigger", NamelistType::BOOLEAN);
  t_nml.addKeyword("ngbr_placement", NamelistType::BOOLEAN);
  t_nml.addKeyword("track_purge", NamelistType::BOOLEAN);
  t_nml.addKeyword("enforce_sanity", NamelistType::BOOLEAN);
  t_nml.addKeyword("max_reports", NamelistType::INTEGER);
  t_nml.addKeyword("interval_trigger", NamelistType::INTEGER,
                   std::to_string(default_inspection_interval));
  t_nml.addKeyword("force_threshold", NamelistType::REAL,
                   realToString(default_large_force_threshold, 9, 4, NumberFormat::STANDARD_REAL));
  t_nml.addKeyword("speed_threshold", NamelistType::REAL,
                   realToString(default_high_speed_threshold, 9, 4, NumberFormat::STANDARD_REAL));
  t_nml.addKeyword("bond_sane", NamelistType::REAL,
                   realToString(default_bond_max_strain, 9, 4, NumberFormat::STANDARD_REAL));
  t_nml.addKeyword("angle_sane", NamelistType::REAL,
                   realToString(default_angle_max_strain, 9, 4, NumberFormat::STANDARD_REAL));
  t_nml.addKeyword("varname", NamelistType::STRING, std::string(default_debug_variable_prefix));

  // Add help messages
  t_nml.addHelp("forces", "Call for checks on the sum of all GPU- or CPU-calculated forces acting "
                "on all atoms and virtual site particles.");
  t_nml.addHelp("force_trigger", "Request that large forces (defined by the value of "
                "force_threshold) be reported and investigated, triggering all active CPU-based "
                "checks on forces.  Without setting this to TRUE, large forces exceeding the "
                "stated threshold can still be recorded for the final report, as long as the run "
                "finishes with success.");
  t_nml.addHelp("ngbr_placement", "Request that the arrangement of particles in the neighbor list "
                "be checked by equivalent CPU-based routines for constucting the cells.");
  t_nml.addHelp("track_purge", "Request that momentum purges, translational as well as rotational "
                "if appropriate, be tracked in the debugging output.");
  t_nml.addHelp("enforce_sanity", "Demand that system sanity checks occur during the course of a "
                "calculation.  Specific checks and consequences of failing the checks will vary "
                "from program to program, but in general systems that fail will be removed from "
                "the calculation with alerts issued to the user.");
  t_nml.addHelp("max_reports", "Set the maximum number of anomalous incidents which will be "
                "recorded.");
  t_nml.addHelp("interval_trigger", "Set the interval where all active checks, whether on forces "
                "or other aspects of the calculation, will be triggered regardless of whether "
                "anything anomalous seems to have appeared.");
  t_nml.addHelp("force_threshold", "Set the magnitude at which a force will be considered large, "
                "in units of kcal/mol-A.");
  t_nml.addHelp("speed_threshold", "Set the rate of movement at which a particle will be "
                "considered to have high speed, in units of A/fs.");
  t_nml.addHelp("bond_sane", "Set the tolerance for a strained bond to be considered part of a "
                "'sane' structure, in terms of kcal/mol-Angstrom^2.");
  t_nml.addHelp("angle_sane", "Set the tolerance for a strained angle to be considered part of a "
                "'sane' structure, in terms of kcal/mol-radian^2.");
  t_nml.addHelp("varname", "Set the prefix for variables in the output report (the file name of "
                "the report is specified in the &files namelist).  This string will be extended "
                "by descriptive modifiers for the results of different debugging checks.");

  // There is expected to be at most one unique &debug namelist in a given input file.  Seek it
  // out by wrapping back to the beginning of the input file if necessary.
  *start_line = readNamelist(tf, &t_nml, *start_line, wrap, tf.getLineCount(), found);

  return t_nml;
}

} // namespace namelist
} // namespace stormm
