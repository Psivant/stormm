// -*-c++-*-
#ifndef CONFORMER_ANALYSIS_H
#define CONFORMER_ANALYSIS_H

#include <vector>
#include "copyright.h"
#include "../../../src/Accelerator/gpu_details.h"
#include "../../../src/MoleculeFormat/mdlmol.h"
#include "../../../src/Namelists/nml_conformer.h"
#include "../../../src/Namelists/nml_files.h"
#include "../../../src/Namelists/nml_minimize.h"
#include "../../../src/Namelists/nml_report.h"
#include "../../../src/Namelists/user_settings.h"
#include "../../../src/Potential/scorecard.h"
#include "../../../src/Potential/static_exclusionmask.h"
#include "../../../src/Synthesis/phasespace_synthesis.h"
#include "../../../src/Synthesis/synthesis_cache_map.h"
#include "../../../src/Synthesis/synthesis_permutor.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/Structure/clash_detection.h"
#include "../../../src/Structure/local_arrangement.h"

namespace conf_app {
namespace analysis {

using stormm::card::GpuDetails;
using stormm::energy::ScoreCard;
using stormm::energy::StaticExclusionMask;
using stormm::namelist::ConformerControls;
using stormm::namelist::default_minimize_clash_r0;
using stormm::namelist::default_minimize_clash_ratio;
using stormm::namelist::FilesControls;
using stormm::namelist::ReportControls;
using stormm::namelist::UserSettings;
using stormm::structure::ClashReport;
using stormm::structure::MdlMol;
using stormm::synthesis::Condensate;
using stormm::synthesis::PhaseSpaceSynthesis;
using stormm::synthesis::PsSynthesisReader;
using stormm::synthesis::SynthesisCacheMap;
using stormm::synthesis::SynthesisPermutor;
using stormm::synthesis::SystemCache;
using stormm::synthesis::SystemGrouping;
    
/// \brief Compute the angles made by each rotatable bond.  This will compute the angle based on
///        the highest Z-number atoms attached to either end of the bond (not counting the atoms at
///        the other end of the bond).  If there is a tie, the atom with the lower topological
///        index will be preferred.

/// \brief Filter structures to obtain the ones with the best energies separated by some minimum
///        positional RMSD threshold.
///
/// \param poly_ps        Synthesis of coordinates, holding all of the energy-minimized structures
/// \param poly_ps_masks  Exclusions for each unique topology within the synthesis
/// \param sc             Cache of systems used to seed the synthesis
/// \param scmap          Map between the systems cache and the working synthesis
/// \param emin           The minimum energy values and history
/// \param confcon        Conformer namelist user input
  
std::vector<int> filterMinimizedStructures(const PhaseSpaceSynthesis &poly_ps,
                                           const std::vector<StaticExclusionMask> &poly_ps_masks,
                                           const SystemCache &sc, const SynthesisCacheMap &scmap,
                                           const ScoreCard &emin, const ConformerControls &confcon,
                                           const GpuDetails &gpu = null_gpu);

/// \brief Print the best structures for each system, grouped by topology, in accord with user
///        input.
///
/// \param poly_ps       Synthesis of topologies from each
/// \param best_confs    List of the best configurations
/// \param emin          Energy tracking with histories from each minimization
/// \param sc            The cache of all systems read from disk
/// \param scmap         Map between the systems cache and the working synthesis
/// \param confcon       Conformer namelist user input
/// \param repcon        Control data from a &report namelist in the user input, to specify
///                      modifications to SD files
void printResults(const PhaseSpaceSynthesis &poly_ps, const std::vector<int> &best_confs,
                  const ScoreCard &emin, const SystemCache &sc, const SynthesisCacheMap &scmap,
                  const ConformerControls &confcon, const ReportControls &repcon);

/// \brief Print a report for the user to summarize the results and indicate ways to improve the
///        outcomes.
///
/// \param sc           Cache of systems used to seed the calculations
/// \param ui           Contains user-generated contol input
/// \param sandbox_prm  Contains information on the permutations and combinatorial space available
///                     to each system in the systems cache sc
/// \param sandbox_map  Map linking the coordinate synthesis used in the initial round of energy
///                     minimizations to the systems cache created based on the &files namelist
/// \param prelim_emin  Outputs from preliminary energy minimization runs
/// \param best_confs   List of the best configurations
void printReport(const SystemCache &sc, const UserSettings &ui,
                 const SynthesisPermutor &sandbox_prm, const SynthesisCacheMap &sandbox_map,
                 const ScoreCard &prelim_emin, const std::vector<int> &best_confs);
  
} // namespace analysis
} // namespace conf_app

#endif
