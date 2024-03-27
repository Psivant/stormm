// -*-c++-*-
#ifndef STORMM_PRESENT_ENERGY_H
#define STORMM_PRESENT_ENERGY_H

#include <fstream>
#include <string>
#include <vector>
#include "copyright.h"
#include "Namelists/nml_files.h"
#include "Namelists/nml_report.h"
#include "Namelists/user_settings.h"
#include "Potential/energy_enumerators.h"
#include "Potential/scorecard.h"
#include "Synthesis/synthesis_cache_map.h"
#include "Synthesis/synthesis_enumerators.h"
#include "report_table.h"

namespace stormm {
namespace review {

using energy::ScoreCard;
using energy::EnergySample;
using energy::StateVariable;
using namelist::FilesControls;
using namelist::ReportControls;
using namelist::UserSettings;
using synthesis::SynthesisCacheMap;
using synthesis::SystemGrouping;

/// \brief Make an assessment of all relevant topologies to determine whether any of the synthesis
///        of systems carry energy terms, if no list is provided a-priori.  Otherwise, simply
///        return the existing list of state variables.
///
/// \param nrg         Collection of tracked energy state variables for all systems in a synthesis.
///                    In principle, the ScoreCard can track energies in any ordered collection of
///                    systems, but the coordinate synthesis provides the framework in this
///                    library.
/// \param scmap       Map linking the synthesis behind nrg to an underlying collectin of user
///                    input (this is how information in an expansive synthesis can make its way
///                    back to the user in terms that they have outlined)
/// \param quantities  The list of energy state variables of interest (if empty, all possible state
///                    variables will be investigated)
std::vector<StateVariable> assessEnergyDetails(const ScoreCard &nrg,
                                               const SynthesisCacheMap &scmap,
                                               const std::vector<StateVariable> &quantities = {});

/// \brief Make a table of the average energies and standard deviations for all requested energy
///        components found in a ScoreCard object.  If no instances of certain energy terms are
///        detected in any of the systems, they will not be reported unless specifically requested.
///
/// \param nrg           The energy outputs of some calculation or series of calculations
/// \param measure       Indicate whether to take averages at each time point in the series,
///                      averages for the final time point only, or averages of averages over each
///                      system's time series.
/// \param varname       Base name of the variable in which energy outputs will be stored
/// \param scmap         Map linking the synthesis behind nrg to an underlying collectin of user
///                      input (this is how information in an expansive synthesis can make its way
///                      back to the user in terms that they have outlined)
/// \param quantities    List of the energy components to report.  An empty vector implies printing
///                      of all quantities so long as they are present in one or more systems.
///                      This input can be used to limit the scope of the output or standardize
///                      outputs from multiple runs to simplify posterior comparisons.
std::vector<ReportTable> tabulateAverageEnergy(const ScoreCard &nrg, EnergySample measure,
                                               const std::string &varname,
                                               const SynthesisCacheMap &scmap,
                                               const std::vector<StateVariable> &quantities = {},
                                               const ReportControls &repcon = ReportControls());

/// \brief Make a table of system designations so that subsequent tables can take indices rather
///        han complete system names when a synthesis is divided into various groupings.
///
///
/// \param scmap         Map linking the synthesis associated with an energy tracking object to
///                      user input.  The user input likely originates in a &files namelist.
/// \param organization  The means of grouping individual systems within the synthesis
/// \param shortcut_key  Pointer to a pre-existing string ready to take an explanation of any path
///                      shortcuts that are found in this function
/// \param repcon        User-specified &report namelist control information
ReportTable groupSystemDesignations(const SynthesisCacheMap &scmap,
                                    const SystemGrouping organization, std::string *shortcut_key,
                                    const ReportControls &repcon);

/// \brief Create a table of one specific energy quantity, averaged over each group of systems.
///        This is called repeatedly by tabulateGroupedEnergy() below until all requested aspects
///        of the energy have been described.  For each group of systems, the mean and standard
///        deviation of the energy quantity of interest are placed in consecutive columns.
///
/// \param nrg             The energy tracking object, containing final states and histories for
///                        each system of a synthesis
/// \param measure         The manner in which to extract energies from the nrg object--time
///                        series, final values, time averages
/// \param system_indices  The indices of systems within the synthesis (and thus the nrg parameter)
///                        comprised by each group
/// \param group_bounds    Bounds array for system_indices
/// \param group_count     The number of relevant groups of systems
/// \param quantity        The energy quantity to tabulate
/// \param varname         Variable name to use in the table (this will be appended with the exact
///                        energy term to differentiate it from other tables)
/// \param repcon          User-specified &report namelist control information
/// \param report_mask     Mask of qualifying systems.  Only those whose indices in this array are
///                        marked TRUE will be reported.  This array must track the system order
///                        of system_indices.
ReportTable groupEnergyAccumulation(const ScoreCard &nrg, const EnergySample measure,
                                    const int* system_indices, const int* group_bounds,
                                    const int group_count, const StateVariable quantity,
                                    const std::string &varname, const ReportControls &repcon);

/// \brief Make a table of the averages and standard deviations of energies found in a ScoreCard
///        object, grouped according to user-specified input details.  This function produces two
///        tables, the first defining indices seen in column headings for the second, where the
///        actual energies are to be found.  The first table holds STRING-type data and will thus
///        be entirely protected behind comment characters.  The second table holds REAL-type data.
///        Descriptions of parameters for this function follow from tabulateAverageEnergy() above,
///        with the addition of:
///
/// \param organization  The means of grouping individual systems within the synthesis behind the
///                      ScoreCard object nrg
/// \param shortcut_key  A key listing path shortcuts used to condense certain topology file names.
///                      Assembled and returned.
/// \param scmap         Map linking the synthesis associated with the energy tracking object to
///                      user input.  The user input likely originates in a &files namelist.
/// \param repcon        User-specified &report namelist control information
/// \param report_mask   A mask of the systems that are to be reported.  If provided as a vector
///                      of booleans, TRUE indicates that the system is valid and should be
///                      reported and all systems must be covered, in the order they appear
///                      according to a direct reading of the list in the systems cache (no
///                      alterations with respect to the group label or original topology).  If
///                      provided as a vector of integers, the system indices of valid systems
///                      should be included.
std::vector<ReportTable> tabulateGroupedEnergy(const ScoreCard &nrg, SystemGrouping organization,
                                               EnergySample measure, const std::string &varname,
                                               std::string *shortcut_key,
                                               const SynthesisCacheMap &scmap,
                                               const std::vector<StateVariable> &quantities = {},
                                               const ReportControls &repcon = ReportControls());

/// \brief Prepare a table assigning numerical identifiers to all systems so that subsequent tables
///        can list a quantity associated with one identifier per column.
///
/// \param scmap         Map linking the synthesis associated with the energy tracking object to
///                      user input.  The user input likely originates in a &files namelist.
/// \param shortcut_key  Pre-existing string to hold the caption indicating any shortcut file paths
///                      found for the table.  Modified and returned.
/// \param repcon        User-specified directives on how to format the report file
ReportTable allSystemDesignations(std::string *shortcut_key, const SynthesisCacheMap &scmap,
                                  const ReportControls &repcon = ReportControls());

/// \brief Make a table of all energies for each system found in a ScoreCard object.  Like the
///        function tabulateGroupedEnergy() above, this function produces two tables, the first
///        defining indices seen in column headings (now, one column for each system in the
///        synthesis) of the second table where the actual energies are to be found.  Descriptions
///        of parameters for this function follow from tabulateAverageEnergy() and
///        tabulateGroupedEnergy() above, with the addition of:
///
/// \param post_script  Pre-existing string (cannot be entered as a null pointer) to hold output
///                     file script commands for processing the tables created by this function.
///                     This string is useful in certain situations when the output will be
///                     interpreted by another program.
/// \{
std::vector<ReportTable> tabulateFullEnergy(const ScoreCard &nrg, EnergySample measure,
                                            const std::string &varname, std::string *shortcut_key,
                                            std::string *post_script,
                                            const SynthesisCacheMap &scmap,
                                            const std::vector<StateVariable> &quantities,
                                            const ReportControls &repcon,
                                            const std::vector<bool> &report_mask);

std::vector<ReportTable> tabulateFullEnergy(const ScoreCard &nrg, EnergySample measure,
                                            const std::string &varname, std::string *shortcut_key,
                                            std::string *post_script,
                                            const SynthesisCacheMap &scmap,
                                            const std::vector<StateVariable> &quantities,
                                            const ReportControls &repcon,
                                            const std::vector<int> &report_mask);

std::vector<ReportTable> tabulateFullEnergy(const ScoreCard &nrg, EnergySample measure,
                                            const std::string &varname, std::string *shortcut_key,
                                            std::string *post_script,
                                            const SynthesisCacheMap &scmap,
                                            const std::vector<StateVariable> &quantities = {},
                                            const ReportControls &repcon = ReportControls());
/// \}

/// \brief Make a table of all energies for each system found in a ScoreCard object.  The output is
///        a series of tables for each outlier discovered, with a separate table (STRING-typed) to
///        list each of the systems and their origins.  Descriptions of parameters for this
///        function follow from tabulateGroupedEnergy() and tabulateFullEnergy() above.
std::vector<ReportTable> tabulateOutlierEnergy(const ScoreCard &nrg, EnergySample measure,
                                               std::string *shortcut_key, std::string *post_script,
                                               const SynthesisCacheMap &scmap,
                                               SystemGrouping organization,
                                               const std::vector<StateVariable> &quantities = {},
                                               const ReportControls &repcon = ReportControls());

/// \brief Print ouptut diagnostics report sections based on user preferences.
///
/// \param nrg    Energy tracking data from a set of molecular mechanics calculations
/// \param scmap  Map to connect systems in the synthesis and associated energy tracking to systems
///               in the user-supplied systems cache
/// \param ui     Collection of user input data
void createDiagnosticReport(const ScoreCard &nrg, const SynthesisCacheMap &scmap,
                            const UserSettings &ui);

} // namespace review
} // namespace stormm

#endif

