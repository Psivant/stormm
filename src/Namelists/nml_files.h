// -*-c++-*-
#ifndef STORMM_NML_FILES_H
#define STORMM_NML_FILES_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/textfile.h"
#include "Trajectory/trajectory_enumerators.h"
#include "input.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using constants::ExceptionResponse;
using constants::ModificationPolicy;
using parse::WrapTextSearch;
using trajectory::CoordinateFileKind;
using trajectory::TrajectoryFusion;

/// \brief Default file names and extensions.  The default input file name varies according to the
///        particular application and is therefore not defined within these libraries.
/// \{
constexpr bool default_filecon_read_all_free = false;
constexpr char default_filecon_topology_name[] = "prmtop";
constexpr char default_filecon_coordinate_name[] = "inpcrd";
constexpr char default_filecon_report_name[] = "md.out";
constexpr char default_filecon_trajectory_name[] = "md.crd";
constexpr char default_filecon_checkpoint_name[] = "md.rst";
constexpr char default_filecon_warnings_name[] = "warn.out";
constexpr char default_filecon_errors_name[] = "err.out";
constexpr char default_filecon_result_fusion[] = "AUTO";
constexpr char default_filecon_sdf_mod_policy[] = "NO";
constexpr char default_filecon_sdf_notification[] = "WARN";
constexpr CoordinateFileKind default_filecon_inpcrd_type = CoordinateFileKind::UNKNOWN;
constexpr CoordinateFileKind default_filecon_outcrd_type = CoordinateFileKind::AMBER_CRD;
constexpr CoordinateFileKind default_filecon_chkcrd_type = CoordinateFileKind::AMBER_ASCII_RST;
constexpr char default_filecon_inpcrd_type_name[] = "AMBER_INPCRD";
constexpr char default_filecon_outcrd_type_name[] = "AMBER_CRD";
constexpr char default_filecon_chkcrd_type_name[] = "AMBER_ASCII_RST";
/// \}

/// \brief Object to encapsulate a system, a coupled set of coordinates and a single topology.
///        This can be used to enforce conf.stormm to read certain coordinates in the context of a
///        particular topology, even when other viable topologies might be available.
class MoleculeSystem {
public:

  /// \brief The constructor can make a blank system or automatically fill values.
  ///
  /// Overloaded:
  ///   - Create a blank object
  ///   - Create an object with a complete list of files and read settings
  ///
  /// \param topology_file_in    Topology file for the system (its type will be detected later)
  /// \param coordinate_file_in  Coordinate file for the system (input)
  /// \param trajectory_file_in  Trajectory file for the system (output)
  /// \param checkpoint_file_in  Restart file to write for the system (output)
  /// \param label_in            System label, i.e. "MyBigProtein"
  /// \param frame_start_in      Starting frame to begin reading input coordinates, if the input
  ///                            coordinate file is a trajectory
  /// \param frame_end_in        Frame at which to stop reading coordinates, if the input
  ///                            coordinate file is a trajectory
  /// \param replica_count_in    Number of replicas of this system to produce in the resulting
  ///                            SystemCache object
  /// \param coordinate_kind_in  File type of the input coordinates file
  /// \param trajectory_kind_in  File type of the resulting trajectory file
  /// \param checkpoint_kind_in  File type of the checkpoint file
  /// \{
  MoleculeSystem();
  MoleculeSystem(const std::string &topology_file_in, const std::string &coordinate_file_in,
                 const std::string &trajectory_file_in, const std::string &checkpoint_file_in,
                 const std::string &label_in, int frame_start_in, int frame_end_in,
                 int replica_count_in, CoordinateFileKind coordinate_kind_in,
                 CoordinateFileKind trajectory_kind_in, CoordinateFileKind checkpoint_kind_in);
  /// \}

  /// \brief This object, containing no const members or pointers to repair, can take the default
  ///        copy and move constructors, plus copy and move assignment operators.
  /// \{
  MoleculeSystem(const MoleculeSystem &original) = default;
  MoleculeSystem(MoleculeSystem &&original) = default;
  MoleculeSystem& operator=(const MoleculeSystem &original) = default;
  MoleculeSystem& operator=(MoleculeSystem &&original) = default;
  /// \}

  /// \brief Get the name of the topology file in this system.
  const std::string& getTopologyFileName() const;

  /// \brief Get the name of the input coordinates file.
  const std::string& getInputCoordinateFileName() const;

  /// \brief Get the name of the trajectory file.
  const std::string& getTrajectoryFileName() const;

  /// \brief Get the name of the checkpoint file to write for this system.
  const std::string& getCheckpointFileName() const;

  /// \brief Get the assigned label for this system.
  const std::string& getLabel() const;

  /// \brief Get the starting frame that this system will read for initial coordinate states.
  int getStartingFrame() const;

  /// \brief Get the last frame that this system will read for initial coordinate states.
  int getFinalFrame() const;

  /// \brief Get the total number of frames, thus initial states, that this system expects.
  int getTotalFrames() const;

  /// \brief Get the replica count for this system, the number of copies of each initial
  ///        coordinate frame that it will spawn.
  int getReplicaCount() const;

  /// \brief Get the type of input coordinate file.
  CoordinateFileKind getInputCoordinateFileKind() const;

  /// \brief Get the type of trajectory file.
  CoordinateFileKind getTrajectoryFileKind() const;

  /// \brief Get the type of checkpoint file.
  CoordinateFileKind getCheckpointFileKind() const;

  /// \brief Set the topology file name.  This is useful if pre-allocating an array of
  ///        MoleculeSystems and then filling it up later.
  void setTopologyFileName(const std::string &file_name);

  /// \brief Set the input coordinates file name.
  ///
  /// \param file_name  The chosen file name
  void setInputCoordinateFileName(const std::string &file_name);

  /// \brief Set the trajectory file name.
  ///
  /// \param file_name  The chosen file name
  void setTrajectoryFileName(const std::string &file_name);

  /// \brief Set the checkpoint file name.
  ///
  /// \param file_name  The chosen file name
  void setCheckpointFileName(const std::string &file_name);

  /// \brief Set the starting frame.  This can be necessary if the input coordinates file does not
  ///        contain the requested numbers of frames.
  ///
  /// \param frame_number  The selected frame number
  void setStartingFrame(int frame_number);

  /// \brief Set the final frame.  This can be necessary if the input coordinates file does not
  ///        contain the requested numbers of frames.
  ///
  /// \param frame_number  The selected frame number
  void setFinalFrame(int frame_number);

  /// \brief Set the number of replicas.
  ///
  /// \param count  The number of system replicas to spawn
  void setReplicaCount(int count);

  /// \brief Set the input coordinates file type.
  ///
  /// Overloaded:
  ///   - Get the file type from a string (will be validated)
  ///   - Get the file type from a direct enumeration
  ///
  /// \param kind  The type of input coordinates file to expect (or that was detected)
  /// \{
  void setInputCoordinateFileKind(const std::string &kind);
  void setInputCoordinateFileKind(CoordinateFileKind kind);
  /// \}

  /// \brief Set the trajectory file type.
  ///
  /// Overloaded:
  ///   - Get the file type from a string (will be validated)
  ///   - Get the file type from a direct enumeration
  ///
  /// \param kind  The type of trajectory file to write
  /// \{
  void setTrajectoryFileKind(const std::string &kind);
  void setTrajectoryFileKind(CoordinateFileKind kind);
  /// \}

  /// \brief Set the checkpoint file type.
  ///
  /// Overloaded:
  ///   - Get the file type from a string (will be validated)
  ///   - Get the file type from a direct enumeration
  ///
  /// \param kind  The type of checkpoint file to write
  /// \{
  void setCheckpointFileKind(const std::string &kind);
  void setCheckpointFileKind(CoordinateFileKind kind);
  /// \}
  
  /// \brief Report whether the topology file named in this system is a valid file.  This validator
  ///        is public so that it can be called by a wholistic validation strategy employed in the
  ///        containing FilesControls object.
  bool validateTopologyFile() const;

  /// \brief Report whether the input coordinate file named in this system is a valid file.
  bool validateInputCoordinateFile() const;
  
private:
  std::string topology_file_name;      ///< Topology file describing the molecular system
  std::string coordinate_file_name;    ///< Coordinate file (may be a whole trajectory)
  std::string coordinate_output_name;  ///< Coordinate output file for this system (if not
                                       ///<   supplied, the name will be inferred based on
                                       ///<   directives in the containing FileControls object)
  std::string checkpoint_name;         ///< Checkpoint file for this system (if not supplied, the
                                       ///<   name will be inferred based on directives in the
                                       ///<   containing FileControls object)
  std::string label;                   ///< Label assigned to this system for future reference
  int frame_start;                     ///< Strating frame of the coordinates file to read
  int frame_end;                       ///< Final frame of the coordinates file to read (the
                                       ///<   program will read frames over the range
                                       ///<   [ start, end ], so [ 0, 0 ] gets the first frame)
  int replica_count;                   ///< Number of times to replicate this system
  CoordinateFileKind coordinate_kind;  ///< Kind of coordinate file to be expected
  CoordinateFileKind trajectory_kind;  ///< Kind of trajectory file to write
  CoordinateFileKind checkpoint_kind;  ///< Kind of checkpoint file to write
};

/// \brief Distill the results of file identification, producing clean lists of free topologies,
///        free coordinate files, and linked topology / coordinate systems.  Record output and
///        trajectory files, if available.
class FilesControls {
public:

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param policy_in         Requested error handling behavior
  /// \param tf                Input file translated into RAM
  /// \param start_line        Line of the input file to begin searching for the &solvent namelist
  /// \param found_nml         Indication that the namelist was found
  /// \param report_name_in    New default name of the report file (for application-specific
  ///                          naming conventions)
  /// \param coord_base_in     New default base name of trajectory or other coordinate files (for
  ///                          application-specific naming conventions)
  /// \param coord_ext_in      New default extension of trajectory or other coordinate files (for
  ///                          application-specific naming conventions)
  /// \param sys_requirements  System requirements specifications.  The -sys keyword will be
  ///                          critical under different circumstances and has no defaults.  Some
  ///                          situations may not require the user to specify all four of the
  ///                          files associated with the -sys keyword, and some situations will
  ///                          require that certain files exist in addition to merely being
  ///                          specified.  This array can hold four critical specifiers, -p, -c,
  ///                          -x, and -r, for "topology", "input coordinates", "trajectory", and
  ///                          "restart / checkpoint" files, respectively, and each can accept the
  ///                          extension "e" to indicate that the file must exist.  Conversely, a
  ///                          "g" extension indicates that, in that particular context, the file
  ///                          subkey is bogus and should trigger an exception consistent with
  ///                          the stated policy (policy_in).  Default requirements are a topology
  ///                          and input coordinates, both of which must exist.
  /// \param wrap              Indicate that the search for a &files namelist should carry on
  ///                          from the beginning of an input file if no such namelist is found
  ///                          starting from the original starting point
  /// \{
  FilesControls(ExceptionResponse policy_in = ExceptionResponse::DIE,
                WrapTextSearch wrap = WrapTextSearch::NO);
  FilesControls(const TextFile &tf, int *start_line, bool *found_nml = nullptr,
                ExceptionResponse policy_in = ExceptionResponse::DIE,
                WrapTextSearch wrap = WrapTextSearch::NO,
                const std::vector<std::string> &alternatives = {},
                const std::vector<std::string> &sys_requirements = {"-pe", "-ce"});
  /// \}

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  /// \{
  FilesControls(const FilesControls &original) = default;
  FilesControls(FilesControls &&original) = default;
  FilesControls& operator=(const FilesControls &original) = default;
  FilesControls& operator=(FilesControls &&original) = default;
  /// \}

  /// \brief Get the structure count, based on the number of free coordinate files as well as the
  ///        number of systems, with frame counts and replicas therein.  If one of the systems
  ///        does not have the requisite number of frames, this is an error and may be trapped.
  ///        The structure count reported here is therefore a maximum that can be expected.
  int getStructureCount() const;

  /// \brief Get the free topology count.
  int getFreeTopologyCount() const;

  /// \brief Get the free coordinate file count.
  int getFreeCoordinatesCount() const;

  /// \brief Get the number of system specifications made with the -sys keyword.
  int getSystemDefinitionCount() const;

  /// \brief Get the indicator of whether to read all free coordinate files' frames
  bool readAllFreeFrames() const;
  
  /// \brief Get the coordinate (trajectory) file output format
  CoordinateFileKind getOutputCoordinateFormat() const;
  
  /// \brief Get the coordinate (checkpoint) file output format
  CoordinateFileKind getCheckpointFormat() const;
  
  /// \brief Get the preferences for fusing output coordinate files.
  TrajectoryFusion getFileFusionProtocol() const;
  
  /// \brief Get one or more free topology names.
  ///
  /// Overloaded:
  ///   - Get the free topology named at one index of the array
  ///   - Get all free topology names in the array
  ///
  /// \param index  The index to query
  /// \{
  std::string getFreeTopologyName(int index) const;
  std::vector<std::string> getFreeTopologyNames() const;
  /// \}

  /// \brief Get one or more free coordinate file names.
  ///
  /// Overloaded:
  ///   - Get the free coordinate file named at one index of the array
  ///   - Get all free coordinate file names in the array
  ///
  /// \param index  The index to query
  /// \{
  std::string getFreeCoordinateName(int index) const;
  std::vector<std::string> getFreeCoordinateNames() const;
  /// \}

  /// \brief Get a molecule system from this object's array.
  ///
  /// \param index  The index to query
  MoleculeSystem getSystem(int index) const;

  /// \brief Get the name of the report file (equivalent to mdout in sander or pmemd)
  std::string getReportFile() const;

  /// \brief Get the name of the input transcript file (there is no sander or pmemd equivalent
  ///        other than the verbose reproduction of the input settings at the beginning of mdout)
  std::string getInputTranscriptFile() const;
  
  /// \brief Get the base name of trajectory files to write
  std::string getTrajectoryFileName() const;

  /// \brief Get the base name of (coordinate) checkpoint files to write
  std::string getCheckpointFileName() const;

  /// \brief Get the name of the file containing warnings printed by the program
  std::string getWarningFileName() const;

  /// \brief Get the policy on modifying .sdf file outputs to conform to the Biovia standard
  ModificationPolicy getSdfModificationPolicy() const;

  /// \brief Get an indication of whether to alert the user when correcting minor mistakes in .sdf
  ///        files.
  ExceptionResponse getSdfNotifications() const;

  /// \brief Get the original namelist emulator object as a transcript of the user input.
  const NamelistEmulator& getTranscript() const;
  
  /// \brief Set whether to read all frames from each free trajectory (true), or just the first
  ///        (false).
  ///
  /// \param active  Indicator of whether to have all free frames read (true), or not (false)
  void setAllFreeFrameReading(const bool active);
  
  /// \brief Set the coordinate (trajectory) file output format
  ///
  /// Overloaded:
  /// - Take a string value specifying the coordinate file kind
  /// - Take the literal, enumerated type
  ///
  /// \param traj_kind  The selected trajectory format
  /// \{
  void setOutputCoordinateFormat(const std::string &traj_kind);
  void setOutputCoordinateFormat(CoordinateFileKind traj_kind);
  /// \}
  
  /// \brief Set the checkpoint (coordinate) file output format
  ///
  /// Overloaded:
  /// - Take a string value specifying the coordinate file kind
  /// - Take the literal, enumerated type
  ///
  /// \param chk_kind  The selected checkpoint format
  /// \{
  void setCheckpointFormat(const std::string &chk_kind);
  void setCheckpointFormat(CoordinateFileKind chk_kind);
  /// \}

  /// \brief Add a free topology to the list, after checking for its prior existence in the list.
  ///
  /// \param file_name  Name of the topology file to add
  void addFreeTopologyName(const std::string &file_name);

  /// \brief Remove entries from the free topologies array.
  ///
  /// Overloaded:
  ///   - Remove topologies in a numbered range
  ///   - Remove a topology by name
  ///
  /// \param index    Starting index for system removal, or the one system to remove
  /// \param stretch  Number of systems, including the starting index, to remove from the list
  /// \{
  void removeFreeTopologyName(int index, int stretch = 1);
  void removeFreeTopologyName(const std::string &fname);
  /// \}
  
  /// \brief Add a free coordinate file to the list, after checking for its prior existence.
  ///
  /// \param file_name  Name of the coordinate file to add
  void addFreeCoordinateName(const std::string &file_name);

  /// \brief Remove entries from the free coordinates array
  ///
  /// Overloaded:
  ///   - Remove coordinate files in a numbered range
  ///   - Remove a coordinate file by name
  ///
  /// \param index    Starting index for system removal, or the one system to remove
  /// \param stretch  Number of systems, including the starting index, to remove from the list
  /// \{
  void removeFreeCoordinateName(int index, int stretch = 1);
  void removeFreeCoordinateName(const std::string &fname);
  /// \}

  /// \brief Add a system to the list.
  ///
  /// \param new_mol  The new molecule system to push to the back of the list
  void addSystem(const MoleculeSystem &new_mol);

  /// \brief Remove entries from the systems array.
  ///
  /// \param index    Starting index for system removal, or the one system to remove
  /// \param stretch  Number of systems, including the starting index, to remove from the list
  void removeSystem(int index, int stretch = 1);

  /// \brief Set the report file name.
  ///
  /// \param file_name  New name for the calculation report file
  void setReportFileName(const std::string &file_name);

  /// \brief Set the input transcript file name.
  ///
  /// \param file_name  New name for the input transcript file
  void setInputTranscriptFileName(const std::string &file_name);
  
  /// \brief Set the general (fallback) coordinate output file name.
  ///
  /// \param proto_name  New generic name for unspecified trajectory files
  void setGeneralTrajectoryFileName(const std::string &proto_name);

  /// \brief Set the general (fallback) checkpoint file name.
  ///
  /// \param proto_name  New generic name for unspecified checkpoint files
  void setGeneralCheckpointFileName(const std::string &proto_name);

  /// \brief Set the warning file name.
  ///
  /// \param file_name  New name for the warnings file
  void setWarningFileName(const std::string &file_name);

  /// \brief Set the .sdf output file modification policy.
  ///
  /// \param policy_in  The new policy to follow when encountering minor issues in .sdf file data
  void setSdfModficiationPolicy(ModificationPolicy policy_in);
  
  /// \brief Set whether to emit notifications about .sdf output file modifications.
  ///
  /// \param policy_in  The policy to follow on making alerts
  void setSdfNotifications(ExceptionResponse policy_in);
  
private:

  /// Action to take when receiving bad input
  ExceptionResponse policy;
  
  // Counts of critical data
  int structure_count;          ///< Total number of initial structures, the sum of all free
                                ///<   coordinate files (for which only the first frame will be 
                                ///<   read, if multiple frames are present) and all frames arising
                                ///<   from system definitions (which can include trajectory files
                                ///<   and span multiple frames)
  int free_topology_count;      ///< The number of free topologies, which will be matched to free
                                ///<   coordinate files
  int free_coordinate_count;    ///< The number of free coordinate files, which will be paired to
                                ///<   free trajectories
  int system_count;             ///< The number of system keyword specifications
  bool all_free_frames;         ///< Flag to have all free coordinate files' frames read (if true),
                                ///<   or just the first (if false, default)
  TrajectoryFusion fuse_files;  ///< Indicate whether to fuse output coordinate files (checkpoint
                                ///<   and trajectory files)

  /// Format of the coordinate output files
  CoordinateFileKind coordinate_input_format;
  
  /// Format of the coordinate output files
  CoordinateFileKind coordinate_output_format;

  /// Format of the coordinate checkpoint files
  CoordinateFileKind coordinate_checkpoint_format;

  /// List of free topologies.  These names can indicate specific files, directories containing
  /// files, or regular expressions.  All such files will be evaluated as possible topologies,
  /// and read as individual objects if they appear promising.
  std::vector<std::string> topology_file_names;

  /// List of free coordinate files or trajectories.  These names can indicate specific files,
  /// directories containing files, or regular expressions.  All such files will be evaluated as
  /// possible sources of coordinates.  If trajectories are supplied, the default behavior will be
  /// to take only the first frame from each trajectory, unless the flag all_free_trajectory_frames
  /// is set to true.
  std::vector<std::string> coordinate_file_names;

  /// List of specific systems.  This will be initialized with user input, if there are any
  /// -sys arguments supplied in the &files namelist.  Afterwards, the list will be extended by
  /// adding any viable topology and coordinate pairs that can be detected from the two lists of
  /// independent topologies and coordinate file names.
  std::vector<MoleculeSystem> systems;

  /// Name of the output file.  This is akin to sander's mdout but much more involved as it
  /// spans all systems.
  std::string report_file;

  /// Name of the file to be used in reporting a transcript of all user input in the context of all
  /// possible input directives, including default values and any omitted options.
  std::string input_transcript_file;
  
  /// Conformation output base name.  Each set of initial coordinates will lead to a trajectory,
  /// whether an evolving series of dynamics snapshots through MD, a set of conformers created by
  /// a series of guided minimizations, or some other protocol.  The coordinate files are
  /// therefore the keys to unique systems--topologies to guide the evolution of the coordinates
  /// are either specified by the user or matched by some automated search for the most reasonable
  /// descriptions.  If the output coordinates for a system are specified, that name will be
  /// honored, but if not then the name below will be used to assemble a reasonable name for the
  /// result based on the name of the initial coordinates file.
  std::string coordinate_output_name;

  /// Checkpoint file base name, i.e. the restart file in Amber's sander or pmemd.  The naming
  /// conventions follow those of coordinate_output_name.
  std::string checkpoint_name;

  /// Warnings output file name
  std::string warning_file_name;

  /// Indication of whether to correct minor errors of input .sdf files, e.g. syntax mistakes in
  /// the names of data items (".sdf properties") when printing output .sdf files
  ModificationPolicy sdf_mod_policy;

  /// Indication of whether to warn the user when minor corrections are made to bring .sdf file
  /// outputs back into alignment with the Biovia standard.
  ExceptionResponse sdf_mod_alert;

  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;
};
  
/// \brief Produce a namelist for specifying basic input and output files, which can take the place
///        of a great deal of command line input in the Amber pmemd and sander programs.
///
/// \param tf                 Input text file to scan immediately after the namelist is created
/// \param start_line         Line at which to begin scanning the input file for the namelist
///                           (this will wrap back to the beginning of the file in search of a
///                           unique &files namelist)
/// \param found              Indication that the namelist was found in the input file
/// \param sys_keyword_reqs   Requirements for the subkeys of the -sys keyword (constructed based
///                           on input to the FilesControls object that calls this function)
/// \param policy             Reaction to exceptions encountered during namelist reading
/// \param wrap               Indicate that the search for a &files namelist should carry on
///                           from the beginning of an input file if no such namelist is found
///                           starting from the original starting point
/// \param crd_input_format   The default type for input coordinate files describing each system
/// \param crd_output_format  The default type for trajectory files written based on each system
/// \param crd_chkpt_format   The default type for restart files written based on each system
NamelistEmulator
filesInput(const TextFile &tf, int *start_line, bool *found,
           const std::vector<KeyRequirement> &sys_keyword_reqs,
           ExceptionResponse policy = ExceptionResponse::DIE,
           WrapTextSearch wrap = WrapTextSearch::NO,
           CoordinateFileKind crd_input_format = default_filecon_inpcrd_type,
           CoordinateFileKind crd_output_format = default_filecon_outcrd_type,
           CoordinateFileKind crd_checkpoint_format = default_filecon_chkcrd_type);

} // namespace namelist
} // namespace stormm

#endif
