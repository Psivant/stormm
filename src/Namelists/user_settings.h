// -*-c++-*-
#ifndef STORMM_USER_SETTINGS_H
#define STORMM_USER_SETTINGS_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "FileManagement/file_enumerators.h"
#include "Namelists/command_line_parser.h"
#include "Namelists/nml_analysis.h"
#include "Namelists/nml_debug.h"
#include "Namelists/nml_dynamics.h"
#include "Namelists/nml_emulate.h"
#include "Namelists/nml_ffmorph.h"
#include "Namelists/nml_files.h"
#include "Namelists/nml_minimize.h"
#include "Namelists/nml_pppm.h"
#include "Namelists/nml_precision.h"
#include "Namelists/nml_random.h"
#include "Namelists/nml_receptor.h"
#include "Namelists/nml_remd.h"
#include "Namelists/nml_report.h"
#include "Namelists/nml_restraint.h"
#include "Namelists/nml_solvent.h"
#include "Parsing/textfile.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/trajectory_enumerators.h"
#include "nml_conformer.h"

namespace stormm {
namespace namelist {

using constants::ExceptionResponse;
using diskutil::PrintSituation;
using parse::TextFile;
using topology::AtomGraph;
using topology::ImplicitSolventModel;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFileKind;

/// \brief Default settings for general file output
/// \{
constexpr PrintSituation default_file_writing_directive = PrintSituation::OPEN_NEW;
/// \}
  
/// \brief Default input settings for various STORMM app input files
/// \{
constexpr char default_conformer_input_file[] = "cgen.in";
constexpr char default_ffrefine_input_file[] = "ffld.in";
/// \}
  
/// \brief Object to hold general user input data, including file names or regular expressions for
///        topology and coordinate files, energy minimization settings, analysis protocols, and
///        output controls.
struct UserSettings {

  /// \brief The constructor requires an input file, similar to mdin for the Amber sander program.
  ///
  /// \param clip      User input obtained from the command line
  /// \param sys_reqs  A list of state descriptors for system files.  This will be passed down to
  ///                  the &files namelist.  To specify -p, -c, -x, or -r in an element of this
  ///                  array means that the topology, input coordinates file, output trajectory
  ///                  file, or checkpoint file must be named.  To append "e" after any of these
  ///                  strings implies that the named file must also exist.  To append "g" after
  ///                  any of these strings implies that to specify such a file name is a bogus
  ///                  input.
  UserSettings(const CommandLineParser &clip,
               const std::vector<std::string> &sys_reqs = { "-pe", "-ce" });

  /// \brief With no const members and Standard Template Library objects comprising the only
  ///        complexity beyond scalar data types, the default copy and move constructors as well
  ///        as copy and move assignment operators can take effect.
  ///
  /// \param original  A pre-existing object to copy or move
  /// \param other     Object in the left hand side of an assignment statement
  /// \{
  UserSettings(const UserSettings &original) = default;
  UserSettings(UserSettings &&original) = default;
  UserSettings& operator=(const UserSettings &original) = default;
  UserSettings& operator=(UserSettings &&original) = default;
  /// \}
  
  /// \brief Get the policy (for passing to other operations that may trigger an exception).
  ExceptionResponse getExceptionBehavior() const;

  /// \brief Get the name of the input file.
  const std::string& getInputFileName() const;

  /// \brief Detect whether a &files namelist was present
  bool getFilesPresence() const;

  /// \brief Detect whether a &debug namelist was present
  bool getDebugPresence() const;

  /// \brief Detect whether a &minimize namelist was present
  bool getMinimizePresence() const;

  /// \brief Detect whether a &solvent namelist was present
  bool getSolventPresence() const;
  
  /// \brief Detect whether a &random namelist was present
  bool getRandomPresence() const;
  
  /// \brief Detect whether a &precision namelist was present
  bool getPrecisionPresence() const;
  
  /// \brief Detect whether a &conformer namelist was present
  bool getConformerPresence() const;

  /// \brief Detect whether a &pppm namelist was present
  bool getPPPMPresence() const;

  /// \brief Detect whether a &dynamics namelist was present
  bool getDynamicsPresence() const;

  /// \brief Detect whether a &remd namelist was present
  bool getRemdPresence() const;
  
  /// \brief Detect whether an &analysis namelist was present
  bool getAnalysisPresence() const;
  
  /// \brief Detect whether an &ffmorph namelist was present
  bool getFFMorphPresence() const;

  /// \brief Detect whether an &emulator namelist was present
  bool getEmulatorPresence() const;
  
  /// \brief Detect whether a &report namelist was present
  bool getReportPresence() const;
  
  /// \brief Get the block of information associated with the &files namelist.
  const FilesControls& getFilesNamelistInfo() const;

  /// \brief Get the block of information associated with the &debug namelist.
  const DebugControls& getDebugNamelistInfo() const;
  
  /// \brief Get the block of information associated with the &minimize namelist.
  const MinimizeControls& getMinimizeNamelistInfo() const;
  
  /// \brief Get the block of information associated with the &solvent namelist.
  const SolventControls& getSolventNamelistInfo() const;
  
  /// \brief Get the block of information associated with the &random namelist.
  const RandomControls& getRandomNamelistInfo() const;

  /// \brief Get the block of information associated with the &precision namelist.
  const PrecisionControls& getPrecisionNamelistInfo() const;

  /// \brief Get the block of information associated with the &conformer namelist.
  const ConformerControls& getConformerNamelistInfo() const;

  /// \brief Get the block of information associated with the &receptor namelist.
  const ReceptorControls& getReceptorNamelistInfo() const;

  /// \brief Get the block of information associated with the &pppm namelist.
  const PPPMControls& getPPPMNamelistInfo() const;
  
  /// \brief Get the block of information associated with the &conformer namelist.
  const DynamicsControls& getDynamicsNamelistInfo() const;

  /// \brief Get the block of information associated with the &conformer namelist.
  const RemdControls& getRemdNamelistInfo() const;

  /// \brief Get force field hyperparameter optimization controls through the &ffmorph namelist.
  const FFMorphControls& getFFMorphNamelistInfo() const;

  /// \brief Get non-bonded emulation potential controls through the &emulator namelist.
  const EmulatorControls& getEmulatorNamelistInfo() const;
  
  /// \brief Get the user-specified diagnostics report features.
  const ReportControls& getReportNamelistInfo() const;

  /// \brief Get the number of &restraint namelist control blocks found in the input.
  int getRestraintNamelistCount() const;
  
  /// \brief Get information from the &restraint namelist control block(s).
  ///
  /// Overloaded:
  ///   - Get the complete list of objects derived from &restraint namelist control blocks
  ///   - Get the interpreted form of a specific restraint namelist control block (the index will
  ///     be checked for validity)
  ///
  /// \param index  The index of the namelist control block to pull from the available list
  ///               (checked for validity)
  /// \{
  const std::vector<RestraintControls>& getRestraintNamelistInfo() const;
  const RestraintControls& getRestraintNamelistInfo(int index) const;
  /// \}

  /// \brief Get the number of &analysis namelist control blocks found in the input.
  int getAnalysisNamelistCount() const;
  
  /// \brief Get information from the &analysis namelist control blocks.  Overloading and
  ///        descriptions of input variables follow from getRestraintNamelistInfo(), above.
  /// \{
  const std::vector<AnalysisControls>& getAnalysisNamelistInfo() const;
  const AnalysisControls& getAnalysisNamelistInfo(int index) const;
  /// \}
    
  /// \brief Produce the file overwriting policy.
  PrintSituation getPrintingPolicy() const;

private:

  ExceptionResponse policy;     ///< Action in the event of bad input
  PrintSituation print_policy;  ///< Policy to take with regard to general output files
  bool has_files_nml;           ///< Indicate the presence of a &files namelist in the input
  bool has_debug_nml;           ///< Indicate the presence of a &debug namelist in the input
  bool has_minimize_nml;        ///< Indicate the presence of a &minimize namelist in the input
  bool has_solvent_nml;         ///< Indicate the presence of a &solvent namelist in the input
  bool has_random_nml;          ///< Indicate the presence of a &random namelist in the input
  bool has_precision_nml;       ///< Indicate the presence of a &precision namelist in the input
  bool has_conformer_nml;       ///< Indicate the presence of a &conformer namelist in the input
  bool has_receptor_nml;        ///< Indicate the presence of a &receptor namelist in the input
  bool has_pppm_nml;            ///< Indicate the presence of a &pppm namelist in the input
  bool has_dynamics_nml;        ///< Indicate the presence of a &dynamics namelist in the input
  bool has_remd_nml;            ///< Indicate the presence of a &remd namelist in the input
  bool has_analysis_nml;        ///< Indicate the presence of an &analysis namelist in the input
  bool has_ffmorph_nml;         ///< Indicate the presence of an &ffmorph namelist in the input
  bool has_emulator_nml;        ///< Indicate the presence of an &emulator namelist in the input
  bool has_report_nml;          ///< Indicate the presence of a &report namelist in the input
  int restraint_nml_count;      ///< Number of &restraint namelists found in the input
  int analysis_nml_count;       ///< Number of &analysis namelists found in the input
  
  /// Name of the original input file
  std::string input_file;

  // Control parameters: these structs each encapsulate their own namelist from the input file.
  FilesControls file_io_input;      ///< All input and output file names, save for the command file
  DebugControls debug_io_input;     ///< Debugging user inputs, used to orchestrate functions which
                                    ///<   run outside the usual flow of various workflows to check
                                    ///<   critical results
  MinimizeControls line_min_input;  ///< Line minimization directives
  SolventControls solvent_input;    ///< Implicit solvent specifications
  RandomControls prng_input;        ///< Random number generator specifications
  PrecisionControls prec_input;     ///< Precision model specifications, including accumulation bit
                                    ///<   settings as well as calculation floating point types
  ConformerControls conf_input;     ///< Conformer generation instructions
  ReceptorControls receptor_input;  ///< Grid-based rigid receptor representation instructions
  PPPMControls pppm_input;          ///< Instructions specific to particle-particle / particle-mesh
                                    ///<   evaluations of long-ranged potentials
  DynamicsControls dyna_input;      ///< Molecular dynamics instructions
  RemdControls remd_input;          ///< Replica exchange molecular dynamics instructions
  FFMorphControls ffmod_input;      ///< Force field modification instructions
  EmulatorControls emul_input;      ///< Non-bonded potential emulation and fitting
  ReportControls diagnostic_input;  ///< Diagnostics report file details and layout
  
  /// There can be many restraint controls sections in the input.  This vector holds them all.
  std::vector<RestraintControls> rstr_inputs;

  /// There can be multiple analysis controls sections in the input.  This vector holds them all.
  std::vector<AnalysisControls> analysis_inputs;
};

} // namespace namelist
} // namespace stormm

#endif
