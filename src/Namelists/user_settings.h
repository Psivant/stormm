// -*-c++-*-
#ifndef STORMM_USER_SETTINGS_H
#define STORMM_USER_SETTINGS_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "FileManagement/file_enumerators.h"
#include "Namelists/nml_dynamics.h"
#include "Namelists/nml_ffmorph.h"
#include "Namelists/nml_files.h"
#include "Namelists/nml_minimize.h"
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

/// \brief Enumerate the various apps using STORMM libraries, and which use this UserSettings
///        object as a way to collect common namelists and other command line input.  This
///        enumerator makes it possible to tailor the contents of UserSettings to a specific
///        application, and to tune the behavior in response to particular inputs.
enum class AppName {
  CONFORMER,  ///< The STORMM conformer generator
  DYNAMICS,   ///< The STORMM molecular dynamics program
  FFREFINE    ///< The STORMM force field refinement tool
};
  
/// \brief Object to hold general user input data, including file names or regular expressions for
///        topology and coordinate files, energy minimization settings, analysis protocols, and
///        output controls.
struct UserSettings {

  /// \brief The constructor requires an input file, similar to mdin for the Amber sander program.
  ///
  /// \param argc  Number of command-line variables
  /// \param argv  List of command line argument strings
  UserSettings(int argc, const char* argv[], AppName prog_set);

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

  /// \brief Get the list of command-line arguments
  const std::vector<std::string>& getCommandLineArguments() const;
  
  /// \brief Detect whether a &files namelist was present
  bool getFilesPresence() const;

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

  /// \brief Detect whether a &dynamics namelist was present
  bool getDynamicsPresence() const;

  /// \brief Detect whether a &rmed namelist was present
  bool getRemdPresence() const;
  
  /// \brief Detect whether an &ffmorph namelist was present
  bool getFFMorphPresence() const;
  
  /// \brief Detect whether a &report namelist was present
  bool getReportPresence() const;
  
  /// \brief Get the block of information associated with the &files namelist.
  const FilesControls& getFilesNamelistInfo() const;
  
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

  /// \brief Get the block of information associated with the &conformer namelist.
  const DynamicsControls& getDynamicsNamelistInfo() const;

  /// \brief Get the block of information associated with the &conformer namelist.
  const RemdControls& getRemdNamelistInfo() const;

  /// \brief Get force field hyperparameter optimization controls through the &ffmorph namelist.
  const FFMorphControls& getFFMorphNamelistInfo() const;

  /// \brief Get the user-specified diagnostics report features.
  const ReportControls& getReportNamelistInfo() const;
  
  /// \brief Get a const reference to the vector of &restraint namelist objects.
  const std::vector<RestraintControls>& getRestraintNamelistInfo() const;

  /// \brief Get one block of information associated with a particular &restraint namelist.
  const RestraintControls& getRestraintNamelistInfo(int index) const;
    
  /// \brief Produce the file overwriting policy.
  PrintSituation getPrintingPolicy() const;

private:

  ExceptionResponse policy;     ///< Action in the event of bad input
  PrintSituation print_policy;  ///< Policy to take with regard to general output files
  bool has_files_nml;           ///< Indicate the presence of a &files namelist in the input
  bool has_minimize_nml;        ///< Indicate the presence of a &minimize namelist in the input
  bool has_solvent_nml;         ///< Indicate the presence of a &solvent namelist in the input
  bool has_random_nml;          ///< Indicate the presence of a &random namelist in the input
  bool has_precision_nml;       ///< Indicate the presence of a &precision namelist in the input
  bool has_conformer_nml;       ///< Indicate the presence of a &conformer namelist in the input
  bool has_receptor_nml;        ///< Indicate the presence of a &receptor namelist in the input
  bool has_dynamics_nml;        ///< Indicate the presence of a &dynamics namelist in the input
  bool has_remd_nml;            ///< Indicate the presence of a &remd namelist in the input
  bool has_ffmorph_nml;         ///< Indicate the presence of an &ffmorph namelist in the input
  bool has_report_nml;          ///< Indicate the presence of a &report namelist in the input
  int restraint_nml_count;      ///< Number of &restraint namelists found in the input
  
  /// Name of the original input file
  std::string input_file;

  /// A record of the command-line arguments
  std::vector<std::string> command_line_args;
  
  // Control parameters: these structs each encapsulate their own namelist from the input file.
  FilesControls file_io_input;      ///< All input and output file names, save for the command file
  MinimizeControls line_min_input;  ///< Line minimization directives
  SolventControls solvent_input;    ///< Implicit solvent specifications
  RandomControls prng_input;        ///< Random number generator specifications
  PrecisionControls prec_input;     ///< Precision model specifications, including accumulation bit
                                    ///<   settings as well as calculation floating point types
  ConformerControls conf_input;     ///< Conformer generation instructions
  ReceptorControls receptor_input;  ///< Grid-based rigid receptor representation instructions
  DynamicsControls dyna_input;      ///< Molecular dynamics instructions
  RemdControls remd_input;          ///< Replica exchange molecular dynamics instructions
  FFMorphControls ffmod_input;      ///< Force field modification instructions
  ReportControls diagnostic_input;  ///< Diagnostics report file details and layout
  
  /// There can be many restraint controls sections in the input.  This vector holds them all.
  std::vector<RestraintControls> rstr_inputs;
};

} // namespace namelist
} // namespace stormm

#endif
